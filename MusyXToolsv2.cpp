// MusyXToolsv2.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <regex>
#include <math.h>
#include <byteswap.h>

extern "C"{
#include <vgmstream/vgmstream.h>
}

#include "libsf2/sf2.hpp"

/** For Reference see:
 * http://hcs64.com/files/DSPADPCM.us.pdf
 * http://www.metroid2002.com/retromodding/wiki/AGSC_(File_Format)
 */

using namespace std;

/* Structures for Sound Archives and their sub-sections */

typedef   signed char      s8;
typedef unsigned char      u8;
typedef   signed short     s16;
typedef unsigned short     u16;
typedef   signed int       s32;
typedef unsigned int       u32;


static inline u32 ReadBE(FILE *f, s32 b)
{
    u32 v = 0;
    for (s32 i = b - 8; i >= 0; i -= 8) v |= fgetc(f) << i;
    return v;
}

static inline u32 ReadLE(FILE *f, u32 b)
{
    u32 v = 0;
    for (u32 i = 0; i<b; i += 8) v |= fgetc(f) << i;
    return v;
}

static inline u32 ZeroPadding(FILE*f, s32 b)
{
    u32 ret = ReadBE(f, b);
    if(ret!=0)
    {
        printf ("Warning at offset %ld: Expected zero padding, got %d\n",ftell(f), ret);
    }

    return ret;
}

typedef struct _dsp_sample
{
    u32 id;
    u32 offset;
    u32 size;
    u32 sampOffset;
    u32 coeffOffset;
    u8 baseNote;
    u8 loopFlag;
    u32 sampleRate;

    /* according to  http://www.metroid2002.com/retromodding/wiki/AGSC_(File_Format)
     * 0: DSP compressed sample
     * 1: DSP compressed drum-sample
     * 2: PCM
     * 3: VADPCM compressed (N64 legacy)
     */
    u8 sampleFormat;

    u32 sampleCount; // no. of raw, i.e. uncompressed pcm frames
    s16* pcm=nullptr;

    u32 loopStart;
    u32 loopLength;
    u32 infoOffset; // adpcm decoder
    u16 adpcmCoeff[16];
    // DSP decoder initial state
    u16 bytesPerFrame;
    u8 ps;           // predictor/scale, actually 2 bytes, higher byte is zero though
    u16 yn1;          // sample history
    u16 yn2;          // sample history

    // DSP decoder loop context
    u8 lps;          // predictor/scale for loop context, actually 2 bytes, higher byte is zero though
    u16 lyn1;         // sample history (n-1) for loop context
    u16 lyn2;         // sample history (n-2) for loop context

    ~_dsp_sample()
    {
        delete this->pcm;
        this->pcm = nullptr;
    }
} dsp;

typedef struct
{
    u32 size;
    u16 id;	//
    u16 sampleID;
    u8 rootKey;
    u8 loopFlag;
    u16 adsrID;
    int adsrIndex;
} macro;

typedef struct
{
    u16 id;
    bool exists;
    u16 sampleID;
    u8 baseNote;
    u8 loopFlag;
    u8 startNote;
    u8 endNote;
    bool adsr;
    double attack;
    double decay;
    double sustain;
    double release;
    s8 transpose;
    u8 volume;
    u8 pan;
    u8 sourr;
    u8 prioOfs;
} noteRegion; // noteRegion

typedef struct
{
    bool exists;
    u32 size;
    u16 id;	//
    u16 dummyA;
    u32 noteCount;
    vector<noteRegion> notes;
} instrument;

typedef struct
{
    bool exists;
    u32 size;
    u16 id;	//
    u16 dummyA;
    u32 noteCount;
    vector<noteRegion> notes;
} layer;

typedef struct
{
    u32 size;
    u16 id;
    u8 attack;
    u8 attackDecimal;
    s32 attackTimecents; // attack time in timecents
    u32 attackMs; // attack time in milliseconds
    u8 decay;
    u8 decayDecimal;
    s32 decayTimecents;
    u32 decayMs;
    u8 sustain;
    u8 sustainDecimal;
    double sustaindB;
    u8 release;
    u8 releaseDecimal;
    u32 releaseMs;
    double releaseTime;
} table;

double LogB(double n, double b)
{
    // log(n)/log(2) is log2.
    return log(n) / log(b);
}
double getPan(double pan)
{
    double realPan = (pan - 64) * 7.936507937;
    if (realPan < 0)
        return realPan + 65536;
    else
        return realPan;
}

double getVolume(double volume)
{
    return 200 * abs(LogB(pow((volume / 127), 2), 10));
}

float getSustain(float sustain)
{
    if (sustain == 0)
        return 900;
    else
        return 200 * abs(LogB(((float)sustain / 100), 10));
}

float timeToTimecents(float time)
{
    float timeCent = floor(1200 * LogB(time, 2));
    if (timeCent > -12000 && timeCent < 0)
        return timeCent + 65536;
    else if (timeCent < -12000)
        return -12000;
    else
        return timeCent;
}

// General MIDI instrument names
static const char *const general_MIDI_instr_names[128] =
{
    "Acoustic Grand Piano", "Bright Acoustic Piano", "Electric Grand Piano", "Honky-tonk Piano", "Rhodes Piano", "Chorused Piano",
    "Harpsichord", "Clavinet", "Celesta", "Glockenspiel", "Music Box", "Vibraphone", "Marimba", "Xylophone", "Tubular Bells", "Dulcimer",
    "Drawbar Organ", "Percussive Organ", "Rock Organ", "Church Organ", "Reed Organ", "Accordion", "Harmonica", "Tango Accordion",
    "Acoustic Guitar (nylon)", "Acoustic Guitar (steel)", "Electric Guitar (jazz)", "Electric Guitar (clean)", "Electric Guitar (muted)",
    "Overdriven Guitar", "Distortion Guitar", "Guitar Harmonics", "Acoustic Bass", "Electric Bass (finger)", "Electric Bass (pick)",
    "Fretless Bass", "Slap Bass 1", "Slap Bass 2", "Synth Bass 1", "Synth Bass 2", "Violin", "Viola", "Cello", "Contrabass",
    "Tremelo Strings", "Pizzicato Strings", "Orchestral Harp", "Timpani", "String Ensemble 1", "String Ensemble 2", "SynthStrings 1",
    "SynthStrings 2", "Choir Aahs", "Voice Oohs", "Synth Voice", "Orchestra Hit", "Trumpet", "Trombone", "Tuba", "Muted Trumpet",
    "French Horn", "Brass Section", "Synth Brass 1", "Synth Brass 2", "Soprano Sax", "Alto Sax", "Tenor Sax", "Baritone Sax",
    "Oboe", "English Horn", "Bassoon", "Clarinet", "Piccolo", "Flute", "Recorder", "Pan Flute", "Bottle Blow", "Shakuhachi", "Whistle",
    "Ocarina", "Lead 1 (square)", "Lead 2 (sawtooth)", "Lead 3 (calliope lead)", "Lead 4 (chiff lead)", "Lead 5 (charang)",
    "Lead 6 (voice)", "Lead 7 (fifths)", "Lead 8 (bass + lead)", "Pad 1 (new age)", "Pad 2 (warm)", "Pad 3 (polysynth)", "Pad 4 (choir)",
    "Pad 5 (bowed)", "Pad 6 (metallic)", "Pad 7 (halo)", "Pad 8 (sweep)", "FX 1 (rain)", "FX 2 (soundtrack)", "FX 3 (crystal)",
    "FX 4 (atmosphere)", "FX 5 (brightness)", "FX 6 (goblins)", "FX 7 (echoes)", "FX 8 (sci-fi)", "Sitar", "Banjo", "Shamisen", "Koto",
    "Kalimba", "Bagpipe", "Fiddle", "Shanai", "Tinkle Bell", "Agogo", "Steel Drums", "Woodblock", "Taiko Drum", "Melodic Tom",
    "Synth Drum", "Reverse Cymbal", "Guitar Fret Noise", "Breath Noise", "Seashore", "Bird Tweet", "Telephone Ring", "Helicopter",
    "Applause", "Gunshot"
};

void extract_data(FILE* src, FILE* dst, int size)
{
#define MMMAX 4096
    uint8_t data[MMMAX];
    int read_max = MMMAX;

    int left = size;
    while(left)
    {
        if(read_max > left)
            read_max = left;

        int actual_read = fread(data, 1, read_max, src);

        if(actual_read == 0)
            break; // EOF

        fwrite(data, actual_read, 1, dst);

        left -= read_max;
    }
}

int samples_to_nibbles(int samples)
{
    int whole_frames = samples / 14;
    int remainder = samples % 14;

    if(remainder > 0)
        return (whole_frames * 16) + remainder + 2;
    else
        return whole_frames * 16;
}

int samples_to_bytes(int samples)
{
    int nibbles = samples_to_nibbles(samples);
    return (nibbles / 2) + (nibbles % 2);
}

int main(int argc, const char* argv[])
{
    int actualInstCount = 0, actualDrumCount = 0, curInstrument = 0;
    instrument instruments[128];
    instrument drums[128];
    vector<instrument> layers;
    u8 tempChar;
    u16 tempID;
    u32 tempSize, tempOffset, nextOffset;

    if (argc < 5)
    {
        printf("Usage:\n%s inst.proj inst.pool inst.sdir inst.samp", argv[0]);
        return 1;
    }

    // Now getting sample info
    FILE *sdir = fopen(argv[3], "rb");
    vector<dsp> dsps;
    fseek(sdir, 0, SEEK_END);
    u32 sdirSize = ftell(sdir);
    fseek(sdir, 0, SEEK_SET);
    //int dspCount = (sdirSize - 4) / (0x20 + 0x28);
    //dsps.resize(dspCount);
    for (unsigned int i = 0; i<=0xFFFF; i++)
    {
        u16 id = ReadBE(sdir, 16);
        u16 padding = ZeroPadding(sdir, 16);

        if((u32)(id<<16 | padding) == 0xFFFFFFFFU)
        {
            // terminator symbol
            break;
        }

        dsp sample;
        sample.id = id;
        printf("Reading sample %X\n", id);

        sample.sampOffset = ReadBE(sdir, 32);

        ZeroPadding(sdir, 32);

        sample.baseNote = ReadBE(sdir, 8);
        ZeroPadding(sdir, 8); // this might be part of the samplerate, since it influences the pitch
        sample.sampleRate = ReadBE(sdir, 16);

        sample.sampleFormat = ReadBE(sdir, 8);
        sample.sampleCount = ReadBE(sdir, 24);

        sample.loopStart = ReadBE(sdir, 32);

        sample.loopLength = ReadBE(sdir, 32);

        if (sample.loopLength > 0)
            sample.loopFlag = 1;
        else
            sample.loopFlag = 0;
        sample.infoOffset = ReadBE(sdir, 32);


        dsps.push_back(sample);
    }

    for (unsigned int i = 0; i < dsps.size(); i++)
    {
        fseek(sdir, dsps[i].infoOffset, SEEK_SET);

        // they are big endian, leave them big, since we dont need them

        fread(&dsps[i].bytesPerFrame, 2,1, sdir); // is this of any use??
        fread(&dsps[i].ps, 1,1, sdir);
        fread(&dsps[i].lps, 1,1, sdir);

        // these four bytes must be the loop history, because they are always zero for non-looped samples, i.e. as it should be according to docs
        // no sure about the order though (n-1 first? n-2 first?)
        fread(&dsps[i].lyn2, 2,1, sdir);
        fread(&dsps[i].lyn1, 2,1, sdir);

        for(int j=0; j<16; j++)
            fread(&dsps[i].adpcmCoeff[j], 2, 1, sdir);
    }

    fclose(sdir);

    // open samp file and write out dsp files
    FILE* samp = fopen(argv[4], "rb");

    for (unsigned int i = 0; i < dsps.size(); i++)
    {
        fseek(samp, dsps[i].sampOffset, SEEK_SET);

        char dsp_path[50];
        sprintf(dsp_path, "%05d (0x%04X).dsp", dsps[i].id, dsps[i].id);


        FILE* dsp = fopen(dsp_path, "wb");
// write standard dsp header

        if (dsps[i].sampleCount > 0xDFFFFFFF)// 0xDFFFFFFF samples = 0xFFFFFFFF nibbles
        {
            printf("skipping dsp %d since it has too many samples", dsps[i].sampleCount);
            continue;
        }

        bool loop_flag;
        int loop_start, loop_end;
        if(dsps[i].loopFlag && dsps[i].loopStart + dsps[i].loopLength <= dsps[i].sampleCount)
        {
            loop_flag = 1;
            loop_start = samples_to_nibbles(dsps[i].loopStart);
            loop_end = samples_to_nibbles(dsps[i].loopStart + dsps[i].loopLength) - 1;
        }
        else
        {
            loop_flag = 0;
            loop_start = 2;// # As per the DSPADPCM docs: "If not looping, specify 2, which is the top sample."
            loop_end = 0;
        }

        struct
        {
            // for header generation during decode
            u32 num_samples;      // total number of RAW samples
            u32 num_adpcm_nibbles; // number of ADPCM nibbles (including frame headers)
            u32 sample_rate;      // Sample rate, in Hz
            // DSP addressing and decode context
            u16 loop_flag;    // 1=LOOPED, 0=NOT LOOPED
            u16 format;       // Always 0x0000, for ADPCM
            u32 sa;           // Start offset address for looped samples (zero for non-looped)
            u32 ea;           // End offset address for looped samples
            u32 ca;           // always zero
            u16 coef[16];     // decode coefficients (eight pairs of 16-bit words)
            // DSP decoder initial state
            u16 gain;         // always zero for ADPCM
            u16 ps;           // predictor/scale
            u16 yn1;          // sample history
            u16 yn2;          // sample history
            // DSP decoder loop context
            u16 lps;          // predictor/scale for loop context
            u16 lyn1;         // sample history (n-1) for loop context
            u16 lyn2;         // sample history (n-2) for loop context
            u16 pad[11] = {0};      // reserved
        } sDSPADPCM;

        sDSPADPCM.num_samples = bswap_32(dsps[i].sampleCount);
        sDSPADPCM.num_adpcm_nibbles = bswap_32(samples_to_nibbles(dsps[i].sampleCount));
        sDSPADPCM.sample_rate = bswap_32(dsps[i].sampleRate);

        sDSPADPCM.loop_flag = bswap_16(loop_flag);
        sDSPADPCM.format = bswap_16(0);
        sDSPADPCM.sa = bswap_32(loop_start);
        sDSPADPCM.ea = bswap_32(loop_end);
        sDSPADPCM.ca = bswap_32(0);

        memcpy(sDSPADPCM.coef, dsps[i].adpcmCoeff, sizeof(sDSPADPCM.coef));

        sDSPADPCM.gain = bswap_16(0);

        // store a byte bigendian in u16
        sDSPADPCM.ps = dsps[i].ps << 8;

        sDSPADPCM.yn1 = dsps[i].yn1;
        sDSPADPCM.yn2 = dsps[i].yn2;

        sDSPADPCM.lps = dsps[i].lps << 8;

        sDSPADPCM.lyn1 = dsps[i].lyn1;
        sDSPADPCM.lyn2 = dsps[i].lyn2;

        // write header
        fwrite(&sDSPADPCM, sizeof(sDSPADPCM), 1, dsp);

        int sample_size = samples_to_bytes(dsps[i].sampleCount);
        extract_data(samp, dsp, sample_size);

        fclose(dsp);


        VGMSTREAM * handle = init_vgmstream(dsp_path);
        dsps[i].pcm = new s16[handle->num_samples];
        render_vgmstream(dsps[i].pcm, handle->num_samples, handle);
        close_vgmstream (handle);
    }

    fclose(samp);


    // Reading pool
    FILE *pool = fopen(argv[2], "rb");
    fseek(pool, 0, SEEK_END);
    u32 poolSize = ftell(pool);
    fseek(pool, 0, SEEK_SET);
    u32 macroOffset = ReadBE(pool, 32);
    u32 adsrOffset = ReadBE(pool, 32);
    u32 keymapOffset = ReadBE(pool, 32);
    u32 layerOffset = ReadBE(pool, 32);

    // Checking ADSR first
    int tableCount = 0;
    vector<table> tables;
    if(adsrOffset == 0)
    {
        // seen in SFA (retail) starfoxs.poo
        printf("that pool doesnt contain adsr data\n");
    }
    else
    {
        printf("Checking ADSR tables\n");
        fseek(pool, adsrOffset, SEEK_SET);
        nextOffset = tempOffset = ftell(pool);
        while (ftell(pool) < keymapOffset - 4)
        {
            // size of this table chunk
            tempSize = ReadBE(pool, 32);
            nextOffset += tempSize;

            // ID of this table
            tempID = ReadBE(pool, 16);
            if (tempID != 0xffff)
            {
                tableCount++;
                printf("Table 0x%X: ", tempID);
                tables.resize(tableCount);
                tables[tableCount - 1].size = tempSize;
                tables[tableCount - 1].id = tempID;
                ZeroPadding(pool, 16);

                if (tempSize == 0x1c)
                {
                    printf("contains a DLS ADSR envelope\n");

                    //fseek(pool, 2, SEEK_CUR);

                    tables[tableCount - 1].attackTimecents = ReadLE(pool, 32);
                    //fseek(pool, 2, SEEK_CUR);
                    tables[tableCount - 1].decayTimecents = ReadLE(pool, 32);
                    tables[tableCount - 1].sustaindB = (double)(0x1000 - ReadLE(pool, 16)) * 0.025;
                    tables[tableCount - 1].releaseMs = ReadLE(pool, 16);

                    printf("\tAttack = %d timecents:\n\tDecay = %d timecents:\n\tSustain Level = -%.03f dB:\n\tRelease = %d milliseconds:\n", tables[tableCount - 1].attackTimecents, tables[tableCount - 1].decayTimecents, tables[tableCount - 1].sustaindB, tables[tableCount - 1].releaseMs);

                    // velocity to attack scale
                    // explanation needed what this means
                    ReadLE(pool, 32);

                    // key to decay scale
                    // ???
                    ReadLE(pool, 32);
                }
                else if(tempSize == 0x10)
                {
                    printf("contains a ADSR envelope\n");

                    tables[tableCount - 1].attackMs/*Timecents*/ = ReadLE(pool, 16);
                    tables[tableCount - 1].decayMs/*Timecents*/ = ReadLE(pool, 16);
                    tables[tableCount - 1].sustaindB = (double)(0x1000 - ReadLE(pool, 16)) * 0.025;
                    tables[tableCount - 1].releaseMs = ReadLE(pool, 16);
                    printf("\tAttack = %d milliseconds:\n\tDecay = %d milliseconds:\n\tSustain Level = -%.03f dB:\n\tRelease = %d milliseconds:\n", tables[tableCount - 1].attackMs, tables[tableCount - 1].decayMs, tables[tableCount - 1].sustaindB, tables[tableCount - 1].releaseMs);

                    /*
                    tables[tableCount - 1].attack = ReadBE(pool, 8);
                    tables[tableCount - 1].attackDecimal = ReadBE(pool, 8);
                    tables[tableCount - 1].decay = ReadBE(pool, 8);
                    tables[tableCount - 1].decayDecimal = ReadBE(pool, 8);
                    tables[tableCount - 1].sustain = ReadBE(pool, 8);
                    tables[tableCount - 1].sustainDecimal = ReadBE(pool, 8);
                    tables[tableCount - 1].release = ReadBE(pool, 8);
                    tables[tableCount - 1].releaseDecimal = ReadBE(pool, 8);
                    */
                }
                else
                {
                    printf("contains unknown garbage!\n");
                }
            }
            // seek to the next table chunk
            fseek(pool, nextOffset, SEEK_SET);
        }
    }


    // Should be at macro offset, but let's seek to be sure
    fseek(pool, macroOffset, SEEK_SET);
    nextOffset = tempOffset = ftell(pool);
    vector<macro> macros;
    int macroCount = 0;
    while (ftell(pool) < adsrOffset)
    {
        tempSize = ReadBE(pool, 32);
        nextOffset += tempSize;
        tempID = ReadBE(pool, 16); // SoundMacro ObjectID
        if (tempID != 0xffff)
        {
            macroCount++;
            macros.resize(macroCount);
            macros[macroCount - 1].id = tempID;

            fseek(pool, 2, SEEK_CUR);
            tempOffset = ftell(pool);

            // Looking for sample info (i.e. SoundMacros)
            while (ftell(pool) < nextOffset)
            {
                fseek(pool, 3, SEEK_CUR);
                tempChar = ReadBE(pool, 8);
                if (tempChar == 0x10)
                {
                    fseek(pool, -3, SEEK_CUR);
                    macros[macroCount - 1].sampleID = ReadBE(pool, 16);
                    for (unsigned int i = 0; i < dsps.size(); i++)
                    {
                        if (dsps[i].id == macros[macroCount - 1].sampleID)
                        {
                            macros[macroCount - 1].rootKey = dsps[i].baseNote;
                            macros[macroCount - 1].loopFlag = dsps[i].loopFlag;
                        }
                    }
                    printf("Macro %X uses sample # %X\n", macros[macroCount - 1].id, macros[macroCount - 1].sampleID);
                    fseek(pool, 5, SEEK_CUR);
                }
                else if (tempChar == 0xc)
                {
                    fseek(pool, -3, SEEK_CUR);
                    macros[macroCount - 1].adsrID = ReadBE(pool, 16);
                    for (int i = 0; i < tableCount; i++)
                    {
                        if (tables[i].id == macros[macroCount - 1].adsrID)
                        {
                            macros[macroCount - 1].adsrIndex = i;
                        }
                    }
                    fseek(pool, 5, SEEK_CUR);
                }
                else
                    fseek(pool, 4, SEEK_CUR);
            }
        }
        else
            fseek(pool, nextOffset, SEEK_CUR);

    }

    // Checking layers
    printf("Checking instrument layers\n");
    fseek(pool, layerOffset, SEEK_SET);
    nextOffset = tempOffset = ftell(pool);
    int layerCount = 0;
    while (ftell(pool) < poolSize - 12)
    {
        tempSize = ReadBE(pool, 32);
        if (tempSize != 0xffffffff)
        {
            layerCount++;
            layers.resize(layerCount);
            nextOffset += tempSize;
            tempID = ReadBE(pool, 16);
            fseek(pool, 2, SEEK_CUR);
            layers[layerCount - 1].id = tempID;

            layers[layerCount - 1].noteCount = ReadBE(pool, 32);
            layers[layerCount - 1].notes.resize(layers[layerCount - 1].noteCount);
            printf("Layer %X at 0x%X with %d note regions\n", layers[layerCount - 1].id, ftell(pool), layers[layerCount - 1].noteCount);
            for (unsigned int j = 0; j < layers[layerCount - 1].noteCount; j++)
            {
                tempID = ReadBE(pool, 16);
                if (tempID == 0xffff)
                {
                    fseek(pool, 10, SEEK_CUR);
                    layers[layerCount - 1].notes[j].exists = false;
                }
                else
                {
                    for (int k = 0; k < macroCount; k++)
                    {
                        if (macros[k].id == tempID)
                        {
                            printf("\tNote region %d uses macro %X\n", j, macros[k].id);
                            layers[layerCount - 1].notes[j].sampleID = macros[k].sampleID;
                            layers[layerCount - 1].notes[j].baseNote = macros[k].rootKey;
                            layers[layerCount - 1].notes[j].loopFlag = macros[k].loopFlag;
                            layers[layerCount - 1].notes[j].exists = true;
                            layers[layerCount - 1].notes[j].startNote = ReadBE(pool, 8);
                            layers[layerCount - 1].notes[j].endNote = ReadBE(pool, 8);
                            layers[layerCount - 1].notes[j].transpose = ReadBE(pool, 8);
                            layers[layerCount - 1].notes[j].volume = ReadBE(pool, 8);
                            fseek(pool, 2, SEEK_CUR); // skip priority and surround pan
                            layers[layerCount - 1].notes[j].pan = ReadBE(pool, 8);
                            fseek(pool, 3, SEEK_CUR);
                            if (macros[k].adsrIndex != 0)
                            {
                                layers[layerCount - 1].notes[j].adsr = true;
                                layers[layerCount - 1].notes[j].attack = tables[macros[j].adsrIndex].attackTimecents;
                                layers[layerCount - 1].notes[j].decay = tables[macros[j].adsrIndex].decayTimecents;
                                layers[layerCount - 1].notes[j].sustain = tables[macros[j].adsrIndex].sustaindB * 10;
                                layers[layerCount - 1].notes[j].release = timeToTimecents(tables[macros[j].adsrIndex].releaseTime);
                                /*
                                layers[layerCount - 1].notes[j].attack = (float)(tables[macros[j].adsrIndex].attack / 1e3) + (float)(tables[macros[j].adsrIndex].attackDecimal * 256 / 1e6);
                                layers[layerCount - 1].notes[j].decay = (float)(tables[macros[j].adsrIndex].decay / 1e3) + (float)(tables[macros[j].adsrIndex].decayDecimal * 256 / 1e6);
                                layers[layerCount - 1].notes[j].sustain = (float)(tables[macros[k].adsrIndex].sustain) * 0.0244 + (float)(tables[macros[k].adsrIndex].sustainDecimal * 6.25);
                                layers[layerCount - 1].notes[j].release = (float)(tables[macros[j].adsrIndex].release / 1e3) + (float)(tables[macros[j].adsrIndex].releaseDecimal * 256 / 1e6);
                                */
                            }
                            break;
                        }
                        //						else
                        //							fseek(pool, 10, SEEK_CUR);
                    }
                }
            }
        }
    }

    // Checking drum tables
    printf("Checking keymaps\n");
    fseek(pool, keymapOffset, SEEK_SET);
    int keymapCount = (layerOffset - keymapOffset - 4) / 0x408;
    layerCount += keymapCount;
    layers.resize(layerCount);
    for (int i = layerCount - keymapCount; i < layerCount; i++)
    {
        printf("Reading keymap at 0x%X\n", ftell(pool));
        fseek(pool, 4, SEEK_CUR);
        tempID = ReadBE(pool, 16);
        fseek(pool, 2, SEEK_CUR);
        layers[i].id = tempID;
        layers[i].noteCount = 128;
        layers[i].notes.resize(128);

        // Now to try and map all note regions
        for (int j = 0; j < 128; j++)
        {
            layers[i].notes[j].exists = false;
//				printf("Reading Keymap %d: Note Region %d\n", curInstrument, j);
            tempID = ReadBE(pool, 16);
            if (tempID == 0xffff)
            {
                fseek(pool, 6, SEEK_CUR);
            }
            else if (tempID & 0x8000)  	// Maps to layer, rather than macro
            {
                for (int k = 0; k < layerCount - 1; k++)  	// Reading all layers prior to this one
                {
                    if (layers[k].id == tempID)
                    {
                        printf("Keymap %X note %d uses layer %X\n", layers[i].id, j, layers[k].id);
                        layers[i].noteCount += layers[k].noteCount - 1;
                        layers[i].notes.resize(layers[i].noteCount);
                        layers[i].notes[j] = layers[k].notes[0];	// Only taking the first region for now
                        layers[i].notes[j].transpose = ReadBE(pool, 8);
                        layers[i].notes[j].pan = ReadBE(pool, 8);
                        layers[i].notes[j].exists = true;
                        for (unsigned int n = 1; n < layers[k].noteCount; n++)
                        {
                            layers[i].notes[layers[i].noteCount - n] = layers[k].notes[n];
                            layers[i].notes[layers[i].noteCount - n].transpose = layers[i].notes[j].transpose;
                            layers[i].notes[layers[i].noteCount - n].pan = layers[i].notes[j].pan;
                            layers[i].notes[layers[i].noteCount - n].exists = true;
                            layers[i].notes[layers[i].noteCount - n].startNote = j;
                            layers[i].notes[layers[i].noteCount - n].endNote = j;
                        }
                        fseek(pool, 4, SEEK_CUR);
                        break;
                    }
                }

            }
            else
            {
                for (int k = 0; k < macroCount; k++)
                {
                    if (macros[k].id == tempID)
                    {
                        printf("Keymap %X note %d uses macro %X\n", layers[i].id, j, macros[k].id);
                        layers[i].notes[j].sampleID = macros[k].sampleID;
                        layers[i].notes[j].baseNote = macros[k].rootKey;
                        layers[i].notes[j].loopFlag = macros[k].loopFlag;
                        if (macros[k].adsrIndex != 0)
                        {
                            layers[i].notes[j].adsr = true;
                            layers[i].notes[j].attack = tables[macros[j].adsrIndex].attackTimecents;
                            layers[i].notes[j].decay = tables[macros[j].adsrIndex].decayTimecents;
                            layers[i].notes[j].sustain = tables[macros[j].adsrIndex].sustaindB * 10;
                            layers[i].notes[j].release = timeToTimecents(tables[macros[j].adsrIndex].releaseTime);

                            /*
                            layers[i].notes[j].attack = (float)(tables[macros[j].adsrIndex].attack / 1e3) + (float)(tables[macros[j].adsrIndex].attackDecimal * 256 / 1e6);
                            layers[i].notes[j].decay = (float)(tables[macros[j].adsrIndex].decay / 1e3) + (float)(tables[macros[j].adsrIndex].decayDecimal * 256 / 1e6);
                            layers[i].notes[j].sustain = (float)(tables[macros[k].adsrIndex].sustain) * 0.0244 + (float)(tables[macros[k].adsrIndex].sustainDecimal * 6.25);
                            layers[i].notes[j].release = (float)(tables[macros[j].adsrIndex].release / 1e3) + (float)(tables[macros[j].adsrIndex].releaseDecimal * 256 / 1e6);
                            */
                        }
                        layers[i].notes[j].exists = true;
                        layers[i].notes[j].transpose = ReadBE(pool, 8);
                        layers[i].notes[j].pan = ReadBE(pool, 8);
                        fseek(pool, 4, SEEK_CUR);
                        break;
                    }
                }
            }
            layers[i].notes[j].startNote = j;
            layers[i].notes[j].endNote = j;
        }
    }
    fclose(pool);


    // Getting our instrument info from proj
    FILE *proj = fopen(argv[1], "rb");
    printf("opening .proj\n");
    u32 nextGroupOffset = ReadBE(proj, 32);
    u16 groupID = ReadBE(proj, 16);
    u16 groupType = ReadBE(proj, 16);
    printf("handling groupID %u of type %u\n",groupID, groupType);

    int sizeofEntry = groupType==0 ? 6 : 0x10; // 0 for SongGroup (ie normal / drum page entries of size 6 for use with CSNG), 1 for SFXGroup.

    fseek(proj, 0x1c-8, SEEK_CUR);

    u32 projInstOffset = ReadBE(proj, 32);
    u32 projDrumOffset = ReadBE(proj, 32);
    u32 projMidiSetupTabOffset = ReadBE(proj, 32);
    int instCount = (projDrumOffset - projInstOffset) / sizeofEntry;
    int drumCount = (projMidiSetupTabOffset - projDrumOffset) / sizeofEntry;
    printf("containing %d instruments and %d drum page entries\n",instCount, drumCount);

    fseek(proj, projInstOffset, SEEK_SET);

    for (unsigned short i = 0; i < 128; i++)
    {
        instruments[i].exists = false;
        instruments[i].noteCount = 0;
        instruments[i].notes.resize(1);
        instruments[i].notes[0].sampleID = -1;
        drums[i].exists = false;
        drums[i].noteCount = 0;
        drums[i].notes.resize(1);
        drums[i].notes[0].sampleID = -1;
    }
    for (int i = 0; i < instCount; i++)
    {
        fseek(proj, projInstOffset + i * sizeofEntry, SEEK_SET);
        tempID = ReadBE(proj, 16); // ObjectID of this normal page entry
        if (tempID == 0xffff)
            continue;
        else if (tempID & 0x8000)  	// Normal layer section
        {
            fseek(proj, 2, SEEK_CUR);// skip voice prio and polyphony
            tempChar = ReadBE(proj, 8); // GM midi program number

            for (int j = 0; j < layerCount; j++)
            {
                if (layers[j].id == tempID)
                {
                    instruments[(int)tempChar] = layers[j];
                    break;
                }
            }
            instruments[(int)tempChar].exists = true;
//				printf("Instrument %d exists\n", tempChar);
            printf("%s exists\n", general_MIDI_instr_names[tempChar]);
        }

        else if (tempID & 0x4000)  	// Keymap section
        {
            fseek(proj, 2, SEEK_CUR);// skip voice prio and polyphony
            tempChar = ReadBE(proj, 8);  // GM midi program number


            for (int j = 0; j < layerCount; j++)
            {
                if (layers[j].id == tempID)
                {
                    instruments[(int)tempChar] = layers[j];
                    break;
                }
            }
            instruments[(int)tempChar].exists = true;
//				printf("Instrument %d exists as a keymap\n", tempChar);
            printf("%s exists as a keymap\n", general_MIDI_instr_names[tempChar]);
        }

        else  	// Instrument just has info at macro
        {
            fseek(proj, 2, SEEK_CUR);// skip voice prio and polyphony
            tempChar = ReadBE(proj, 8);  // GM midi program number
            instruments[(int)tempChar].exists = true;

            instruments[(int)tempChar].id = tempID;
            instruments[(int)tempChar].notes.resize(1);
            instruments[(int)tempChar].notes[0].startNote = 0;
            instruments[(int)tempChar].notes[0].endNote = 127;
//				printf("Instrument %d exists as a single macro\n", tempChar);
            printf("%s exists as a single macro\n", general_MIDI_instr_names[tempChar]);
        }
    }

    fseek(proj, projDrumOffset, SEEK_SET);
    for (int i = 0; i < drumCount; i++)
    {
        fseek(proj, projDrumOffset + i * sizeofEntry, SEEK_SET);
        tempID = ReadBE(proj, 16); // ObjectID of this drum page entry
        if (tempID == 0xffff)
        {
            continue;
        }
        else
        {
            fseek(proj, 2, SEEK_CUR); // skip voice prio and polyphony
            tempChar = ReadBE(proj, 8); // GM midi program number

            for (int j = 0; j < layerCount; j++)
            {
                if (layers[j].id == tempID)
                {
                    drums[(int)tempChar] = layers[j];
                    drums[(int)tempChar].exists = true;
                    break;
                }
            }

            if(!drums[(int)tempChar].exists)
            {
                printf("drumkit %d not found in layers, looking through macros\n", (int)tempChar);
                for (unsigned int j = 0; j < macros.size(); j++)
                {
                    if (macros[j].id == tempID)
                    {
                        drums[(int)tempChar].exists = true;

                        drums[(int)tempChar].id = tempID;
                        drums[(int)tempChar].notes.resize(1);
                        drums[(int)tempChar].notes[0].startNote = 0;
                        drums[(int)tempChar].notes[0].endNote = 127;
                        drums[(int)tempChar].notes[0].baseNote = macros[j].rootKey;
                        drums[(int)tempChar].notes[0].sampleID = macros[j].sampleID;
                        drums[(int)tempChar].notes[0].loopFlag = macros[j].loopFlag;
                        printf("drumkit %d exists as a single macro\n", (int)tempChar);
                        break;
                    }
                }
            }
            printf("Drumkit %d ", tempChar);
            if(drums[(int)tempChar].exists)
            {
                printf(" exist\n");
            }
            else
            {
                printf(" does not exist\n");
            }
        }
    }
    fclose(proj);

    printf("Taking care of instruments with only macros\n");
    for (unsigned short i = 0; i < 128; i++)
    {
        if (instruments[i].exists && instruments[i].notes[0].sampleID == (u16)-1)
        {
            printf("Looking for instrument %d macro\n", i);
            for (int j = 0; j < macroCount; j++)
            {
                if (macros[j].id == instruments[i].id)
                {
                    printf("\tMacro %X\n", macros[j].id);
                    instruments[i].notes[0].sampleID = macros[j].sampleID;
                    instruments[i].notes[0].baseNote = macros[j].rootKey;

                    if (macros[j].adsrIndex != 0)
                    {
                        instruments[i].notes[0].adsr = true;
                        instruments[i].notes[0].attack = tables[macros[j].adsrIndex].attackTimecents;
                        instruments[i].notes[0].decay = tables[macros[j].adsrIndex].decayTimecents;
                        instruments[i].notes[0].sustain = tables[macros[j].adsrIndex].sustaindB * 10;
                        instruments[i].notes[0].release = timeToTimecents(tables[macros[j].adsrIndex].releaseTime);
                        /*
                        instruments[i].notes[0].attack = (float)(tables[macros[j].adsrIndex].attack / 1e3) + (float)(tables[macros[j].adsrIndex].attackDecimal * 256 / 1e6);
                        instruments[i].notes[0].decay = (float)(tables[macros[j].adsrIndex].decay / 1e3) + (float)(tables[macros[j].adsrIndex].decayDecimal * 256 / 1e6);
                        instruments[i].notes[0].sustain = (float)(tables[macros[k].adsrIndex].sustain) * 0.0244 + (float)(tables[macros[k].adsrIndex].sustainDecimal * 6.25);
                        instruments[i].notes[0].release = (float)(tables[macros[j].adsrIndex].release / 1e3) + (float)(tables[macros[j].adsrIndex].releaseDecimal * 256 / 1e6);
                        */
                    }

                    printf("\tRoot Key %X\n", macros[j].rootKey);
                    instruments[i].notes[0].exists = true;
                    instruments[i].notes[0].pan = 64;	// Assuming this sample is centered
                    instruments[i].noteCount = 1;

                    // Reserved for ADSR
                }
            }
        }
    }

    ofstream bankTemplate("soundfontBuild.txt");
    stringstream bankTemplateText;
    string bankText;
    bankTemplateText << "[Samples]\n";

    printf("Writing samples\n");
    for (unsigned int i = 0; i < dsps.size(); i++)
    {
        bankTemplateText << "\n    SampleName=" << hex << dsps[i].id << "\n        SampleRate=" << to_string(dsps[i].sampleRate) << "\n        Key=" << to_string(dsps[i].baseNote) << "\n        FineTune=0\n        Type=1\n";
    }
    bankTemplateText << "\n\n[Instruments]\n";


    for (unsigned short i = 0; i < 128; i++)
    {
        if (drums[i].exists && drums[i].noteCount)
        {

            bankTemplateText << "\n    InstrumentName=Drum" << i << "\n";
            for (unsigned int j = 0; j < drums[i].noteCount; j++)
            {

                if (drums[i].notes[j].exists)
                {
                    printf("Printing Drum %d: Note Region %d\n", i, j);

                    bankTemplateText << "\n        Sample=" << hex << drums[i].notes[j].sampleID << "\n";
                    bankTemplateText << "            Z_LowKey=" << to_string(drums[i].notes[j].startNote) << "\n";
                    bankTemplateText << "            Z_HighKey=" << to_string(drums[i].notes[j].endNote) << "\n";
                    bankTemplateText << "            Z_LowVelocity=0\n";
                    bankTemplateText << "            Z_HighVelocity=127\n";
                    bankTemplateText << "            Z_overridingRootKey=" << to_string(drums[i].notes[j].baseNote - drums[i].notes[j].transpose) << "\n"; // + drums[i].notes[j].transpose
                    bankTemplateText << "            Z_initialAttenuation=" << to_string((int)floor(getVolume(drums[i].notes[j].volume))) << "\n";
                    bankTemplateText << "            Z_pan=" << to_string((int)floor(getPan(drums[i].notes[j].pan))) << "\n";

                    if (drums[i].notes[j].adsr)
                    {
                        bankTemplateText << "            Z_attackVolEnv=" << to_string((int)(drums[i].notes[j].attack)) << "\n";
                        bankTemplateText << "            Z_decayVolEnv=" << to_string((int)(drums[i].notes[j].decay)) << "\n";
                        bankTemplateText << "            Z_sustainVolEnv=" << to_string((int)floor((drums[i].notes[j].sustain))) << "\n";
//							bankTemplateText << "            Z_sustainVolEnv=" << to_string((int)floor(drums[i].notes[j].sustain)) << "\n";
                        bankTemplateText << "            Z_releaseVolEnv=" << to_string((int)(drums[i].notes[j].release)) << "\n";
                    }

                    /*
                    bankTemplateText << "            Z_holdVolEnv=" << to_string((int)instruments[j].notes[k].getHold()) << "\n";
                    */
                    bankTemplateText << "            Z_sampleModes=" << dec << (int)(drums[i].notes[j].loopFlag) << "\n";
                }

                else
                    continue;


            }
        }
        else
            continue;
    }

    for (unsigned short i = 0; i < 128; i++)
    {
        if (instruments[i].exists && instruments[i].noteCount)
        {
            printf("Printing Instrument %d:\n", i);

            bankTemplateText << "\n    InstrumentName=" << general_MIDI_instr_names[i] << "\n";
            for (unsigned int j = 0; j < instruments[i].noteCount; j++)
            {

                if (instruments[i].notes[j].exists)
                {
                    printf("\tNote Region %d\n", j);

                    bankTemplateText << "\n        Sample=" << hex << instruments[i].notes[j].sampleID << "\n";
                    bankTemplateText << "            Z_LowKey=" << to_string(instruments[i].notes[j].startNote) << "\n";
                    bankTemplateText << "            Z_HighKey=" << to_string(instruments[i].notes[j].endNote) << "\n";
                    bankTemplateText << "            Z_LowVelocity=0\n";
                    bankTemplateText << "            Z_HighVelocity=127\n";
                    bankTemplateText << "            Z_overridingRootKey=" << to_string(instruments[i].notes[j].baseNote - instruments[i].notes[j].transpose) << "\n"; // + instruments[i].notes[j].transpose
                    bankTemplateText << "            Z_initialAttenuation=" << to_string((int)floor(getVolume(instruments[i].notes[j].volume))) << "\n";
                    bankTemplateText << "            Z_pan=" << to_string((int)floor(getPan(instruments[i].notes[j].pan))) << "\n";

                    if (instruments[i].notes[j].adsr)
                    {
                        bankTemplateText << "            Z_attackVolEnv=" << to_string((int)(instruments[i].notes[j].attack)) << "\n";
                        bankTemplateText << "            Z_decayVolEnv=" << to_string((int)(instruments[i].notes[j].decay)) << "\n";
                        bankTemplateText << "            Z_sustainVolEnv=" << to_string((int)floor((instruments[i].notes[j].sustain))) << "\n";
//							bankTemplateText << "            Z_sustainVolEnv=" << to_string((int)floor(instruments[i].notes[j].sustain)) << "\n";
                        bankTemplateText << "            Z_releaseVolEnv=" << to_string((int)(instruments[i].notes[j].release)) << "\n";
                    }

                    /*
                    bankTemplateText << "            Z_holdVolEnv=" << to_string((int)instruments[j].notes[k].getHold()) << "\n";

                    */
                    bankTemplateText << "            Z_sampleModes=" << dec << (int)(instruments[i].notes[j].loopFlag) << "\n";
                }

                else
                    continue;
            }



        }
        else
            continue;
    }


    bankTemplateText << "\n\n[Presets]\n";

    for (unsigned short i = 0; i < 128; i++)
    {
        if (drums[i].exists && drums[i].noteCount)
        {
            bankTemplateText << "\n    PresetName=Program" << i << "Drum\n        Bank=128\n        Program=" << i << "\n";
            bankTemplateText << "\n        Instrument=Drum" << i << "\n            L_LowKey=0\n            L_HighKey=127\n            L_LowVelocity=0\n            L_HighVelocity=127\n\n";
        }
        else
            continue;
    }

    for (unsigned short i = 0; i < 128; i++)
    {
        if (instruments[i].exists && instruments[i].noteCount)
        {
            bankTemplateText << "\n    PresetName=" << general_MIDI_instr_names[i] << "\n        Bank=0\n        Program=" << i << "\n";
            bankTemplateText << "\n        Instrument=" << general_MIDI_instr_names[i] << "\n            L_LowKey=0\n            L_HighKey=127\n            L_LowVelocity=0\n            L_HighVelocity=127\n\n";
        }
        else
            continue;
    }

    bankTemplateText << "\n[Info]\nVersion=2.1\nEngine=EMU8000 \nName=" << "golf" << "\nROMName=\nROMVersion=0.0\nDate=\nDesigner=\nProduct=\nCopyright=\nEditor=Awave Studio v10.6  \nComments=\n";



    bankText = bankTemplateText.str();
    bankTemplate << bankText;
    bankTemplate.close();


    SF2 sf(22050);
    char buf[50];
    unsigned int bank;
    unsigned int instr;
    unsigned int j;
    unsigned int sf2SampleNum=0, sf2InstrNum=0, sf2BankNum=0;

    for(bank=0; bank<1; bank++)
    {
        instrument (&instrToUse)[128] = instruments;
        for(instr=0; instr<=0x7f; instr++)
        {
            if(((instr / 128) != 0) && ((instr %128)==0))
            {
                // on every multiple of 128 use a new bank, since a bank can only hold 128 instr (from 0 to 0x7f)
                cout << "BankNum: " << sf2BankNum <<endl;
                sf2BankNum++;
            }

            if (instrToUse[instr].exists && instrToUse[instr].noteCount)
            {
                sprintf (buf, "B%02XI%04X", bank, instrToUse[instr].id);
                sf.add_new_instrument(buf);
                    sf2InstrNum++;

                for(j=0;  j < instrToUse[instr].noteCount; j++, sf2SampleNum++)
                {
                    const int& dspIdx = instrToUse[instr].notes[j].sampleID;

                    sprintf(buf, "%05d (0x%04X).dsp", dsps[dspIdx].id, dsps[dspIdx].id);

                    if (dsps[dspIdx].pcm == nullptr)
                    {
                        printf("dsp sample %d has not been decoded", dspIdx);
                        break;
                    }

                    sprintf (buf, "B%02XI%04XS%04X", bank, instrToUse[instr].id, j);
                    sf.add_new_sample(dsps[dspIdx].pcm,
                                      SampleType::SIGNED_16,
                                      buf,
                                      0x00,
                                      dsps[dspIdx].sampleCount,
                                      dsps[dspIdx].loopFlag,
                                      dsps[dspIdx].loopStart,
                                      dsps[dspIdx].loopStart + dsps[dspIdx].loopLength,
                                      dsps[dspIdx].baseNote,
                                      0, //dsps[dspIdx].fineTune, ???
                                      dsps[dspIdx].sampleRate);


                    if (!instrToUse[instr].notes[j].exists)
                    {
                        continue;
                    }

                    sf.add_new_inst_bag();

                    sf.add_new_inst_generator(SFGenerator::sampleModes, dsps[dspIdx].loopFlag); //looping
                    sf.add_new_inst_generator(SFGenerator::sampleID, sf2SampleNum);
                    sf.add_new_inst_generator(SFGenerator::keyRange, instrToUse[instr].notes[j].startNote, instrToUse[instr].notes[j].endNote); // specify key range
                    sf.add_new_inst_generator(SFGenerator::velRange, 0, 127);
                    sf.add_new_inst_generator(SFGenerator::overridingRootKey, instrToUse[instr].notes[j].baseNote - instrToUse[instr].notes[j].transpose);
                    sf.add_new_inst_generator(SFGenerator::initialAttenuation, (uint16_t)floor(getVolume(instrToUse[instr].notes[j].volume)));
                    sf.add_new_inst_generator(SFGenerator::pan, (uint16_t)floor(getPan(instrToUse[instr].notes[j].pan)));

                    // ADSR
//                    sf.add_new_inst_generator(SFGenerator::attackVolEnv, instrToUse[instr].notes[j].attack);
//                    sf.add_new_inst_generator(SFGenerator::decayVolEnv, instrToUse[instr].notes[j].decay);
//                    sf.add_new_inst_generator(SFGenerator::sustainVolEnv, floor((instrToUse[instr].notes[j].sustain)));
//                    sf.add_new_inst_generator(SFGenerator::releaseVolEnv, instrToUse[instr].notes[j].release);

                }// foreach sample

                sprintf (buf, "B%02XI%04X", bank, instrToUse[instr].id);
                sf.add_new_preset(buf, instr%128, bank);
                sf.add_new_preset_bag();
                sf.add_new_preset_generator(SFGenerator::instrument, sf2InstrNum-1);
            }
        }
    }

    cout << "attempting to write sf2 now..." << endl;
    FILE* sfout = fopen("sftest.sf2", "wb");
    if (!sfout)
    {
        cerr<<"error: cannot open soundfont"<<endl;
        return 2;
    }

    sf.write(sfout);

    return 0;
}
