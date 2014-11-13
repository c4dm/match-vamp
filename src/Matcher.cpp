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

#include <iostream>

#include <cstdlib>
#include <cassert>

bool Matcher::silent = true;

//#define DEBUG_MATCHER 1

Matcher::Matcher(Parameters parameters,
                 FeatureExtractor::Parameters feParams,
                 Matcher *p) :
    params(parameters),
    featureExtractor(feParams),
    metric(parameters.distanceNorm)
{
#ifdef DEBUG_MATCHER
    cerr << "Matcher::Matcher(" << params.sampleRate << ", " << p << ")" << endl;
#endif

    otherMatcher = p;	// the first matcher will need this to be set later
    firstPM = (!p);
    frameCount = 0;
    runCount = 0;
    featureSize = featureExtractor.getFeatureSize();
    blockSize = 0;

    blockSize = lrint(params.blockTime / params.hopTime);
#ifdef DEBUG_MATCHER
    cerr << "Matcher: blockSize = " << blockSize << endl;
#endif

    initialised = false;
}

Matcher::Matcher(Parameters parameters, Matcher *p, int featureSize_) :
    params(parameters),
    featureSize(featureSize_),
    featureExtractor(FeatureExtractor::Parameters(params.sampleRate, params.fftSize)), // unused default config
    metric(parameters.distanceNorm)
{
#ifdef DEBUG_MATCHER
    cerr << "Matcher::Matcher(" << params.sampleRate << ", " << p << ", " << featureSize << ")" << endl;
#endif

    otherMatcher = p;	// the first matcher will need this to be set later
    firstPM = (!p);
    frameCount = 0;
    runCount = 0;
    blockSize = 0;

    blockSize = lrint(params.blockTime / params.hopTime);
#ifdef DEBUG_MATCHER
    cerr << "Matcher: blockSize = " << blockSize << endl;
#endif

    initialised = false;
} 

Matcher::~Matcher()
{
#ifdef DEBUG_MATCHER
    cerr << "Matcher(" << this << ")::~Matcher()" << endl;
#endif
}

void
Matcher::init()
{
    if (initialised) return;

    frames = vector<vector<double> >
        (blockSize, vector<double>(featureSize, 0));

    distXSize = blockSize * 2;
    expand();

    frameCount = 0;
    runCount = 0;
    
    initialised = true;
}

void
Matcher::expand()
{
    int distSize = (params.maxRunCount + 1) * blockSize;

    bestPathCost.resize(distXSize, vector<int>(distSize, 0));
    distance.resize(distXSize, vector<unsigned char>(distSize, 0));
    distYSizes.resize(blockSize, distSize);
    first.resize(distXSize, 0);
    last.resize(distXSize, 0);
}

vector<double>
Matcher::consumeFrame(double *reBuffer, double *imBuffer)
{
    if (!initialised) init();

    vector<double> real(reBuffer, reBuffer + params.fftSize/2 + 1);
    vector<double> imag(imBuffer, imBuffer + params.fftSize/2 + 1);
    vector<double> feature = featureExtractor.process(real, imag);
    int frameIndex = frameCount % blockSize;
    frames[frameIndex] = feature;
    calcAdvance();

    return feature;
}

void
Matcher::consumeFeatureVector(std::vector<double> feature)
{
    if (!initialised) init();
    int frameIndex = frameCount % blockSize; 
    frames[frameIndex] = feature;
    calcAdvance();
}

void
Matcher::calcAdvance()
{
    int frameIndex = frameCount % blockSize;

    if (frameCount >= distXSize) {
        expand();
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
        distance[frameCount - blockSize].resize(len, 0);
        for (int i = 0; i < len; ++i) {
            distance[frameCount - blockSize][i] =
                distance[frameCount][i];
        }

        bestPathCost[frameCount] = bestPathCost[frameCount - blockSize];
        bestPathCost[frameCount - blockSize].resize(len, 0);
        for (int i = 0; i < len; ++i) {
            bestPathCost[frameCount - blockSize][i] =
                bestPathCost[frameCount][i];
        }

        distYSizes[frameCount] = distYSizes[frameCount - blockSize];
        distYSizes[frameCount - blockSize] = len;
    }

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

        int dMN = metric.calcDistanceScaled
            (frames[frameIndex],
             otherMatcher->frames[index % blockSize],
             params.distanceScale);
        
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

    frameCount++;
    runCount++;

    otherMatcher->runCount = 0;

    if (overflow && !silent)
        cerr << "WARNING: overflow in distance metric: "
             << "frame " << frameCount << ", val = " << mx << endl;
    
    if (!silent)
        std::cerr << "Frame " << frameCount << ", d = " << (mx-mn) << std::endl;
}

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
            otherMatcher->bestPathCost[j].resize(idx * 2, 0);
            otherMatcher->distance[j].resize(idx * 2, 0);
        }
        otherMatcher->distance[j][idx] = (unsigned char)((dMN & MASK) | dir);
        otherMatcher->bestPathCost[j][idx] =
            (value + (dir==ADVANCE_BOTH? dMN*2: dMN));
    }
} // setValue()

