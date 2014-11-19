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

FeatureExtractor::FeatureExtractor(Parameters parameters) :
    m_params(parameters),
    m_ltAverage(0)
{
    m_featureSize = getFeatureSizeFor(parameters);
    m_prevFrame = vector<double>(m_featureSize, 0.0);

    makeFreqMap();
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
#ifdef DEBUG_MATCHER
        cerr << "makeFreqMap: calling makeChromaFrequencyMap" << endl;
#endif
        makeChromaFrequencyMap();
    } else {
#ifdef DEBUG_MATCHER
        cerr << "makeFreqMap: calling makeStandardFrequencyMap" << endl;
#endif
        makeStandardFrequencyMap();
    }
}

void
FeatureExtractor::makeStandardFrequencyMap()
{
    double binWidth = m_params.sampleRate / m_params.fftSize;
    int crossoverBin = (int)(2 / (pow(2, 1/12.0) - 1));
    int crossoverMidi = lrint(log(crossoverBin*binWidth/440.0)/
                              log(2.0) * 12 + 69);

    // freq = 440 * Math.pow(2, (midi-69)/12.0) / binWidth;
    
    int i = 0;
    while (i <= crossoverBin) {
        m_freqMap[i] = i;
        ++i;
    }

    while (i <= m_params.fftSize/2) {
        double midi = log(i*binWidth/440.0) / log(2.0) * 12 + 69;
        if (midi > 127) midi = 127;
        int target = crossoverBin + lrint(midi) - crossoverMidi;
        if (target >= m_featureSize) target = m_featureSize - 1;
        m_freqMap[i++] = target;
    }
}

void
FeatureExtractor::makeChromaFrequencyMap()
{
    double binWidth = m_params.sampleRate / m_params.fftSize;
    int crossoverBin = (int)(1 / (pow(2, 1/12.0) - 1));
    int i = 0;
    while (i <= crossoverBin) {
        m_freqMap[i++] = 0;
    }
    while (i <= m_params.fftSize/2) {
        double midi = log(i*binWidth/440.0) / log(2.0) * 12 + 69;
        m_freqMap[i++] = (lrint(midi)) % 12 + 1;
    }
}

vector<double>
FeatureExtractor::process(const vector<double> &real, const vector<double> &imag)
{
    vector<double> frame(m_featureSize, 0.0);
    
    double rms = 0;
    for (int i = 0; i <= m_params.fftSize/2; i++) {
        double mag = real[i] * real[i] + imag[i] * imag[i];
        rms += mag;
        frame[m_freqMap[i]] += mag;
    }
    rms = sqrt(rms / (m_params.fftSize/2));

    return postProcess(frame, rms);
}

vector<double>
FeatureExtractor::process(const float *cframe)
{
    vector<double> frame(m_featureSize, 0.0);
    
    double rms = 0;
    for (int i = 0; i <= m_params.fftSize/2; i++) {
        double mag = cframe[i*2] * cframe[i*2] + cframe[i*2+1] * cframe[i*2+1];
        rms += mag;
        frame[m_freqMap[i]] += mag;
    }
    rms = sqrt(rms / (m_params.fftSize/2));

    return postProcess(frame, rms);
}

vector<double>
FeatureExtractor::postProcess(const vector<double> &frame, double rms)
{
    vector<double> feature(m_featureSize, 0.0);

    double totalEnergy = 0;
    if (m_params.useSpectralDifference) {
        for (int i = 0; i < m_featureSize; i++) {
            totalEnergy += frame[i];
            if (frame[i] > m_prevFrame[i]) {
                feature[i] = frame[i] - m_prevFrame[i];
            } else {
                feature[i] = 0;
            }
        }
    } else {
        for (int i = 0; i < m_featureSize; i++) {
            feature[i] = frame[i];
            totalEnergy += feature[i];
        }
    }

    if (m_ltAverage == 0) {
	m_ltAverage = totalEnergy;
    } else {
	double decay = m_params.decay;
        m_ltAverage = m_ltAverage * decay + totalEnergy * (1.0 - decay);
    }

    if (rms <= m_params.silenceThreshold) {
        for (int i = 0; i < m_featureSize; i++) {
            feature[i] = 0;
	}
    } else if (m_params.frameNorm == NormaliseFrameToSum1) {
        for (int i = 0; i < m_featureSize; i++) { 
            feature[i] /= totalEnergy;
	}
    } else if (m_params.frameNorm == NormaliseFrameToLTAverage) {
        for (int i = 0; i < m_featureSize; i++) {
            feature[i] /= m_ltAverage;
	}
    }

    m_prevFrame = frame;
    return feature;
}
    
