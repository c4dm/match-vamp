/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    This file copyright 2007 Simon Dixon, Chris Cannam and QMUL.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "Matcher.h"
#include "Finder.h"

#include <iostream>

#include <cstdlib>

bool Matcher::silent = true;

const double Matcher::decay = 0.99;
const double Matcher::silenceThreshold = 0.0004;
const int Matcher::MAX_RUN_COUNT = 3;

//#define DEBUG_MATCHER 1

Matcher::Matcher(float rate, Matcher *p)
{
#ifdef DEBUG_MATCHER
    std::cerr << "Matcher::Matcher(" << rate << ", " << p << ")" << std::endl;
#endif

    sampleRate = rate;
    otherMatcher = p;	// the first matcher will need this to be set later
    firstPM = (!p);
    matchFileOffset = 0;
    ltAverage = 0;
    frameCount = 0;
    runCount = 0;
    paused = false;
    hopSize = 0;
    fftSize = 0;
    blockSize = 0;
    hopTime = 0.020;	// DEFAULT, overridden with -h  //!!!
    fftTime = 0.04644;	// DEFAULT, overridden with -f
    blockTime = 10.0;	// DEFAULT, overridden with -c
    normalise1 = true;
    normalise2 = false;
    normalise3 = false;
    normalise4 = true;
    useSpectralDifference = true;
    useChromaFrequencyMap = false;
    scale = 90;
    maxFrames = 0;	// stop at EOF

    hopSize = lrint(sampleRate * hopTime);
    fftSize = lrint(pow(2.0, (int)lrint(log(fftTime * sampleRate) / log(2.0))));
    blockSize = lrint(blockTime / hopTime);

    distance = 0;
    bestPathCost = 0;
    distYSizes = 0;
    distXSize = 0;

    initialised = false;

} // default constructor

void
Matcher::setHopSize(int sz)
{
    if (initialised) {
        std::cerr << "Matcher::setHopSize: Can't set after use" << std::endl;
        return;
    }

    hopSize = sz;
    hopTime = float(hopSize) / sampleRate;
    blockTime = blockSize * hopTime;
}

Matcher::~Matcher()
{
#ifdef DEBUG_MATCHER
    std::cerr << "Matcher(" << this << ")::~Matcher()" << std::endl;
#endif

    if (initialised) {
        
        for (int i = 0; i < distXSize; ++i) {
            if (distance[i]) {
                free(distance[i]);
                free(bestPathCost[i]);
            }
        }
        free(distance);
        free(bestPathCost);

        free(first);
        free(last);

        free(distYSizes);
    }
}

void
Matcher::print()
{
    cerr << toString() << endl;
} // print()

string
Matcher::toString()
{
    std::stringstream os;
    os << "Matcher " << this << ": (" << sampleRate
       << "kHz)"
       << "\n\tHop size: " << hopSize 
       << "\n\tFFT size: " << fftSize 
       << "\n\tBlock size: " << blockSize;
    return os.str();
} // toString()

void
Matcher::init()
{
    if (initialised) return;

    initialised = true;

    makeFreqMap(fftSize, sampleRate);

    initVector<double>(prevFrame, freqMapSize);
    initVector<double>(newFrame, freqMapSize);
    initMatrix<double>(frames, blockSize, freqMapSize);
    initVector<double>(totalEnergies, blockSize);

    int distSize = (MAX_RUN_COUNT + 1) * blockSize;

    distXSize = blockSize * 2;

//    std::cerr << "Matcher::init: distXSize = " << distXSize << std::endl;

    distance = (unsigned char **)malloc(distXSize * sizeof(unsigned char *));
    bestPathCost = (int **)malloc(distXSize * sizeof(int *));
    distYSizes = (int *)malloc(distXSize * sizeof(int));

    for (int i = 0; i < blockSize; ++i) {
        distance[i] = (unsigned char *)malloc(distSize * sizeof(unsigned char));
        bestPathCost[i] = (int *)malloc(distSize * sizeof(int));
        distYSizes[i] = distSize;
    }
    for (int i = blockSize; i < distXSize; ++i) {
        distance[i] = 0;
    }
    
    first = (int *)malloc(distXSize * sizeof(int));
    last = (int *)malloc(distXSize * sizeof(int));

    frameCount = 0;
    runCount = 0;
//    frameRMS = 0;
    ltAverage = 0;

    if (!silent) print();
} // init

void
Matcher::makeFreqMap(int fftSize, float sampleRate)
{
    initVector<int>(freqMap, fftSize/2 + 1);
    if (useChromaFrequencyMap)
        makeChromaFrequencyMap(fftSize, sampleRate);
    else
        makeStandardFrequencyMap(fftSize, sampleRate);
} // makeFreqMap()

void
Matcher::makeStandardFrequencyMap(int fftSize, float sampleRate)
{
    double binWidth = sampleRate / fftSize;
    int crossoverBin = (int)(2 / (pow(2, 1/12.0) - 1));
    int crossoverMidi = lrint(log(crossoverBin*binWidth/440.0)/
                              log(2.0) * 12 + 69);
    // freq = 440 * Math.pow(2, (midi-69)/12.0) / binWidth;
    int i = 0;
    while (i <= crossoverBin) {
        freqMap[i] = i;
        ++i;
    }
    while (i <= fftSize/2) {
        double midi = log(i*binWidth/440.0) / log(2.0) * 12 + 69;
        if (midi > 127)
            midi = 127;
        freqMap[i++] = crossoverBin + lrint(midi) - crossoverMidi;
    }
    freqMapSize = freqMap[i-1] + 1;
    if (!silent) {
        cerr << "Standard map size: " << freqMapSize 
             << ";  Crossover at: " << crossoverBin << endl;
//!!!            for (i = 0; i < fftSize / 2; i++)
//                cerr << "freqMap[" << i << "] = " << freqMap[i] << endl;
    }
} // makeStandardFrequencyMap()

void
Matcher::makeChromaFrequencyMap(int fftSize, float sampleRate)
{
    double binWidth = sampleRate / fftSize;
    int crossoverBin = (int)(1 / (pow(2, 1/12.0) - 1));
    // freq = 440 * Math.pow(2, (midi-69)/12.0) / binWidth;
    int i = 0;
    while (i <= crossoverBin)
        freqMap[i++] = 0;
    while (i <= fftSize/2) {
        double midi = log(i*binWidth/440.0) / log(2.0) * 12 + 69;
        freqMap[i++] = (lrint(midi)) % 12 + 1;
    }
    freqMapSize = 13;
    if (!silent) {
        cerr << "Chroma map size: " << freqMapSize 
             << ";  Crossover at: " << crossoverBin << endl;
        for (i = 0; i < fftSize / 2; i++)
            cerr << "freqMap[" << i << "] = " << freqMap[i] << endl;
    }
} // makeChromaFrequencyMap()

vector<double>
Matcher::processFrame(double *reBuffer, double *imBuffer)
{
    if (!initialised) init();

    for (int i = 0; i < (int)newFrame.size(); ++i) {
        newFrame[i] = 0;
    }
    double rms = 0;
    for (int i = 0; i <= fftSize/2; i++) {
        double mag = reBuffer[i] * reBuffer[i] +
                     imBuffer[i] * imBuffer[i];
        rms += mag;
        newFrame[freqMap[i]] += mag;
    }
    rms = sqrt(rms / (fftSize/2));

    int frameIndex = frameCount % blockSize;

    if (frameCount >= distXSize) {
//        std::cerr << "Resizing " << distXSize << " -> " << distXSize * 2 << std::endl;
        distXSize *= 2;
        distance = (unsigned char **)realloc(distance, distXSize * sizeof(unsigned char *));
        bestPathCost = (int **)realloc(bestPathCost, distXSize * sizeof(int *));
        distYSizes = (int *)realloc(distYSizes, distXSize * sizeof(int));
        first = (int *)realloc(first, distXSize * sizeof(int));
        last = (int *)realloc(last, distXSize * sizeof(int));
        
        for (int i = distXSize/2; i < distXSize; ++i) {
            distance[i] = 0;
        }
    }

    if (firstPM && (frameCount >= blockSize)) {

        int len = last[frameCount - blockSize] -
                 first[frameCount - blockSize];

        // We need to copy distance[frameCount-blockSize] to
        // distance[frameCount], and then truncate
        // distance[frameCount-blockSize] to its first len elements.
        // Same for bestPathCost.
/*
        std::cerr << "Matcher(" << this << "): moving " << distYSizes[frameCount - blockSize] << " from " << frameCount - blockSize << " to "
                  << frameCount << ", allocating " << len << " for "
                  << frameCount - blockSize << std::endl;
*/
        distance[frameCount] = distance[frameCount - blockSize];

        distance[frameCount - blockSize] = (unsigned char *)
            malloc(len * sizeof(unsigned char));
        for (int i = 0; i < len; ++i) {
            distance[frameCount - blockSize][i] =
                distance[frameCount][i];
        }

        bestPathCost[frameCount] = bestPathCost[frameCount - blockSize];

        bestPathCost[frameCount - blockSize] = (int *)
            malloc(len * sizeof(int));
        for (int i = 0; i < len; ++i) {
            bestPathCost[frameCount - blockSize][i] =
                bestPathCost[frameCount][i];
        }

        distYSizes[frameCount] = distYSizes[frameCount - blockSize];
        distYSizes[frameCount - blockSize] = len;
    }

    double totalEnergy = 0;
    if (useSpectralDifference) {
        for (int i = 0; i < freqMapSize; i++) {
            totalEnergy += newFrame[i];
            if (newFrame[i] > prevFrame[i]) {
                frames[frameIndex][i] = newFrame[i] - prevFrame[i];
            } else {
                frames[frameIndex][i] = 0;
            }
        }
    } else {
        for (int i = 0; i < freqMapSize; i++) {
            frames[frameIndex][i] = newFrame[i];
            totalEnergy += frames[frameIndex][i];
        }
    }
    totalEnergies[frameIndex] = totalEnergy;

    double decay = frameCount >= 200 ? 0.99:
        (frameCount < 100? 0: (frameCount - 100) / 100.0);

    if (ltAverage == 0)
        ltAverage = totalEnergy;
    else
        ltAverage = ltAverage * decay + totalEnergy * (1.0 - decay);

// System.err.println(Format.d(ltAverage,4) + " " +
//					Format.d(totalEnergy) + " " +
//					Format.d(frameRMS));

//    std::cerr << "ltAverage: " << ltAverage << ", totalEnergy: " << totalEnergy << ", frameRMS: " << rms << std::endl;

    if (rms <= 0.01)  //!!! silenceThreshold)
        for (int i = 0; i < freqMapSize; i++)
            frames[frameIndex][i] = 0;
    else if (normalise1)
        for (int i = 0; i < freqMapSize; i++)
            frames[frameIndex][i] /= totalEnergy;
    else if (normalise3)
        for (int i = 0; i < freqMapSize; i++)
            frames[frameIndex][i] /= ltAverage;

    vector<double> processedFrame = frames[frameIndex];

    int stop = otherMatcher->frameCount;
    int index = stop - blockSize;
    if (index < 0)
        index = 0;
    first[frameCount] = index;
    last[frameCount] = stop;

    bool overflow = false;
    int mn= -1;
    int mx= -1;
    for ( ; index < stop; index++) {
        int dMN = calcDistance(frames[frameIndex],
                               otherMatcher->frames[index % blockSize]);
        if (mx<0)
            mx = mn = dMN;
        else if (dMN > mx)
            mx = dMN;
        else if (dMN < mn)
            mn = dMN;
        if (dMN >= 255) {
            overflow = true;
            dMN = 255;
        }
        if ((frameCount == 0) && (index == 0))    // first element
            setValue(0, 0, 0, 0, dMN);
        else if (frameCount == 0)                 // first row
            setValue(0, index, ADVANCE_OTHER,
                     getValue(0, index-1, true), dMN);
        else if (index == 0)                      // first column
            setValue(frameCount, index, ADVANCE_THIS,
                     getValue(frameCount - 1, 0, true), dMN);
        else if (index == otherMatcher->frameCount - blockSize) {
            // missing value(s) due to cutoff
            //  - no previous value in current row (resp. column)
            //  - no diagonal value if prev. dir. == curr. dirn
            int min2 = getValue(frameCount - 1, index, true);
            //	if ((firstPM && (first[frameCount - 1] == index)) ||
            //			(!firstPM && (last[index-1] < frameCount)))
            if (first[frameCount - 1] == index)
                setValue(frameCount, index, ADVANCE_THIS, min2, dMN);
            else {
                int min1 = getValue(frameCount - 1, index - 1, true);
                if (min1 + dMN <= min2)
                    setValue(frameCount, index, ADVANCE_BOTH, min1,dMN);
                else
                    setValue(frameCount, index, ADVANCE_THIS, min2,dMN);
            }
        } else {
            int min1 = getValue(frameCount, index-1, true);
            int min2 = getValue(frameCount - 1, index, true);
            int min3 = getValue(frameCount - 1, index-1, true);
            if (min1 <= min2) {
                if (min3 + dMN <= min1)
                    setValue(frameCount, index, ADVANCE_BOTH, min3,dMN);
                else
                    setValue(frameCount, index, ADVANCE_OTHER,min1,dMN);
            } else {
                if (min3 + dMN <= min2)
                    setValue(frameCount, index, ADVANCE_BOTH, min3,dMN);
                else
                    setValue(frameCount, index, ADVANCE_THIS, min2,dMN);
            }
        }
        otherMatcher->last[index]++;
    } // loop for row (resp. column)

    vector<double> tmp = prevFrame;
    prevFrame = newFrame;
    newFrame = tmp;

    frameCount++;
    runCount++;

    otherMatcher->runCount = 0;

    if (overflow && !silent)
        cerr << "WARNING: overflow in distance metric: "
             << "frame " << frameCount << ", val = " << mx << endl;

    if (!silent)
        std::cerr << "Frame " << frameCount << ", d = " << (mx-mn) << std::endl;

    if ((frameCount % 100) == 0) {
        if (!silent) {
            cerr << "Progress:" << frameCount << " " << ltAverage << endl;
        }
    }

    return processedFrame;
} // processFrame()

int
Matcher::calcDistance(const vector<double> &f1, const vector<double> &f2)
{
    double d = 0;
    double sum = 0;
    for (int i = 0; i < freqMapSize; i++) {
        d += fabs(f1[i] - f2[i]);
        sum += f1[i] + f2[i];
    }
    // System.err.print("   " + Format.d(d,3));
    if (sum == 0)
        return 0;
    if (normalise2)
        return (int)(scale * d / sum);	// 0 <= d/sum <= 2
    if (!normalise4)
        return (int)(scale * d);

    // note if this were to be restored, it would have to use
    // totalEnergies vector instead of f1[freqMapSize] which used to
    // store the total energy:
    //	double weight = (5 + Math.log(f1[freqMapSize] + f2[freqMapSize]))/10.0;

    double weight = (8 + log(sum)) / 10.0;
    // if (weight < mins) {
    // 	mins = weight;
    //	System.err.println(Format.d(mins,3) + " " + Format.d(maxs));
    // }
    // if (weight > maxs) {
    // 	maxs = weight;
    //	System.err.println(Format.d(mins,3) + " " + Format.d(maxs));
    // }
    if (weight < 0)
        weight = 0;
    else if (weight > 1)
        weight = 1;
    return (int)(scale * d / sum * weight);
} // calcDistance()

int
Matcher::getValue(int i, int j, bool firstAttempt)
{
    if (firstPM)
        return bestPathCost[i][j - first[i]];
    else
        return otherMatcher->bestPathCost[j][i - otherMatcher->first[j]];
} // getValue()

void
Matcher::setValue(int i, int j, int dir, int value, int dMN)
{
    if (firstPM) {
        distance[i][j - first[i]] = (unsigned char)((dMN & MASK) | dir);
        bestPathCost[i][j - first[i]] =
            (value + (dir==ADVANCE_BOTH? dMN*2: dMN));
    } else {
        if (dir == ADVANCE_THIS)
            dir = ADVANCE_OTHER;
        else if (dir == ADVANCE_OTHER)
            dir = ADVANCE_THIS;
        int idx = i - otherMatcher->first[j];
        if (idx == (int)otherMatcher->distYSizes[j]) {
            // This should never happen, but if we allow arbitrary
            // pauses in either direction, and arbitrary lengths at
            // end, it is better than a segmentation fault.
            std::cerr << "Emergency resize: " << idx << " -> " << idx * 2 << std::endl;
            otherMatcher->distYSizes[j] = idx * 2;
            otherMatcher->bestPathCost[j] =
                (int *)realloc(otherMatcher->bestPathCost[j],
                               idx * 2 * sizeof(int));
            otherMatcher->distance[j] = 
                (unsigned char *)realloc(otherMatcher->distance[j],
                                         idx * 2 * sizeof(unsigned char));
        }
        otherMatcher->distance[j][idx] = (unsigned char)((dMN & MASK) | dir);
        otherMatcher->bestPathCost[j][idx] =
            (value + (dir==ADVANCE_BOTH? dMN*2: dMN));
    }
} // setValue()

