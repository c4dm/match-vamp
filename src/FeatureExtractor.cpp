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

#include "FeatureExtractor.h"

#include <iostream>

#include <cstdlib>
#include <cassert>
#include <cmath>

using namespace std;

//#define DEBUG_FEATURE_EXTRACTOR 1

FeatureExtractor::FeatureExtractor(Parameters parameters) :
    m_params(parameters)
{
    m_featureSize = getFeatureSizeFor(parameters);
    makeFreqMap();

#ifdef DEBUG_FEATURE_EXTRACTOR
    cerr << "*** FeatureExtractor: sampleRate = " << parameters.sampleRate
         << ", useChromaFrequencyMap = " << parameters.useChromaFrequencyMap
         << ", fftSize = " << parameters.fftSize << endl;
#endif
}

int
FeatureExtractor::getFeatureSizeFor(Parameters parameters)
{
    if (parameters.useChromaFrequencyMap) {
	return 13;
    } else {
	return 84;
    }
}

void
FeatureExtractor::makeFreqMap()
{
    m_freqMap = vector<int>(m_params.fftSize / 2 + 1, 0);

    if (m_params.useChromaFrequencyMap) {
#ifdef DEBUG_FEATURE_EXTRACTOR
        cerr << "makeFreqMap: calling makeChromaFrequencyMap" << endl;
#endif
        makeChromaFrequencyMap();
    } else {
#ifdef DEBUG_FEATURE_EXTRACTOR
        cerr << "makeFreqMap: calling makeStandardFrequencyMap" << endl;
#endif
        makeStandardFrequencyMap();
    }
}

void
FeatureExtractor::makeStandardFrequencyMap()
{
    double refFreq = m_params.referenceFrequency;
    double binWidth = m_params.sampleRate / m_params.fftSize;
    int crossoverBin = (int)(2 / (pow(2, 1/12.0) - 1));
    int crossoverMidi = lrint(log(crossoverBin * binWidth / refFreq)/
                              log(2.0) * 12 + 69);

#ifdef DEBUG_FEATURE_EXTRACTOR
    cerr << "FeatureExtractor::makeStandardFrequencyMap: refFreq = " << refFreq << endl;
#endif
    
    int i = 0;
    while (i <= crossoverBin) {
        double freq = i * binWidth;
        if (freq < m_params.minFrequency || freq > m_params.maxFrequency) {
            m_freqMap[i++] = -1;
        } else {
            m_freqMap[i] = i;
            i++;
        }
    }

    while (i <= m_params.fftSize/2) {
        double freq = i * binWidth;
        if (freq < m_params.minFrequency || freq > m_params.maxFrequency) {
            m_freqMap[i++] = -1;
        } else {
            double midi = log(freq / refFreq) / log(2.0) * 12 + 69;
            if (midi > 127) midi = 127;
            int target = crossoverBin + lrint(midi) - crossoverMidi;
            if (target >= m_featureSize) target = m_featureSize - 1;
            m_freqMap[i++] = target;
        }
    }

#ifdef DEBUG_FEATURE_EXTRACTOR
    cerr << "FeatureExtractor: crossover bin is " << crossoverBin << " for midi "
         << crossoverMidi << endl;
    cerr << "FeatureExtractor: map is:" << endl;
    for (i = 0; i <= m_params.fftSize/2; ++i) {
        cerr << i << ": " << m_freqMap[i] << ", ";
    }
    cerr << endl;
#endif
}

void
FeatureExtractor::makeChromaFrequencyMap()
{
    double refFreq = m_params.referenceFrequency;
    double binWidth = m_params.sampleRate / m_params.fftSize;
    int crossoverBin = (int)(1 / (pow(2, 1/12.0) - 1));
    int i = 0;
    while (i <= crossoverBin) {
        double freq = i * binWidth;
        if (freq < m_params.minFrequency || freq > m_params.maxFrequency) {
            m_freqMap[i++] = -1;
        } else {
            m_freqMap[i++] = 0;
        }
    }
    while (i <= m_params.fftSize/2) {
        double freq = i * binWidth;
        if (freq < m_params.minFrequency || freq > m_params.maxFrequency) {
            m_freqMap[i++] = -1;
        } else {
            double midi = log(freq / refFreq) / log(2.0) * 12 + 69;
            m_freqMap[i++] = (lrint(midi)) % 12 + 1;
        }
    }
}

vector<double>
FeatureExtractor::process(const vector<double> &real, const vector<double> &imag)
{
    vector<double> frame(m_featureSize, 0.0);
    
    for (int i = 0; i <= m_params.fftSize/2; i++) {
        double mag = real[i] * real[i] + imag[i] * imag[i];
        int index = m_freqMap[i];
        if (index >= 0) {
            frame[index] += mag;
        }
    }

    return frame;
}

vector<double>
FeatureExtractor::process(const float *cframe)
{
    vector<double> frame(m_featureSize, 0.0);
    
    for (int i = 0; i <= m_params.fftSize/2; i++) {
        double mag = cframe[i*2] * cframe[i*2] + cframe[i*2+1] * cframe[i*2+1];
        int index = m_freqMap[i];
        if (index >= 0) {
            frame[index] += mag;
        }
    }

    return frame;
}

