/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright (c) 2007-2015 Simon Dixon, Chris Cannam, and Queen Mary
    University of London, Copyright (c) 2014-2015 Tido GmbH.
    
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
    // Our handling of the referenceFrequency parameter depends on the
    // frequency map in use.

    // With the chroma frequency map, we use referenceFrequency to set
    // up the chroma bin frequencies when constructing the map, and
    // then just follow the map (without having to refer to
    // referenceFrequency again) when we get the frequency-domain
    // audio.

    // With the standard frequency map, using referenceFrequency to
    // set up the map doesn't work so well -- it only really affects
    // the crossover frequency, and much of the useful information is
    // below that frequency. What we do instead is to ignore the
    // referenceFrequency when creating the map -- setting it up for
    // 440Hz -- and then use it to scale the individual
    // frequency-domain audio frames before applying the map to them.
    
    double refFreq = 440.; // See above -- *not* the parameter!
    double binWidth = double(m_params.sampleRate) / m_params.fftSize;
    int crossoverBin = int(2 / (pow(2, 1/12.0) - 1));
    int crossoverMidi = int(log(crossoverBin * binWidth / refFreq)/
                            log(2.0) * 12 + 69 + 0.5);

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
            int target = crossoverBin + int(midi + 0.5) - crossoverMidi;
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
    double binWidth = double(m_params.sampleRate) / m_params.fftSize;
    int crossoverBin = int(1 / (pow(2, 1/12.0) - 1));
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
            m_freqMap[i++] = (int(midi + 0.5)) % 12 + 1;
        }
    }
}

feature_t
FeatureExtractor::process(const vector<double> &real, const vector<double> &imag)
{
    vector<float> mags(m_params.fftSize/2 + 1, 0.0);

    for (int i = 0; i <= m_params.fftSize/2; i++) {
        mags[i] = float(real[i] * real[i] + imag[i] * imag[i]);
    }

    return processMags(mags);
}

feature_t
FeatureExtractor::process(const vector<float> &real, const vector<float> &imag)
{
    vector<float> mags(m_params.fftSize/2 + 1, 0.0);

    for (int i = 0; i <= m_params.fftSize/2; i++) {
        mags[i] = real[i] * real[i] + imag[i] * imag[i];
    }

    return processMags(mags);
}

feature_t
FeatureExtractor::process(const float *real, const float *imag)
{
    vector<float> mags(m_params.fftSize/2 + 1, 0.0);

    for (int i = 0; i <= m_params.fftSize/2; i++) {
        mags[i] = real[i] * real[i] + imag[i] * imag[i];
    }

    return processMags(mags);
}

feature_t
FeatureExtractor::process(const float *cframe)
{
    vector<float> mags(m_params.fftSize/2 + 1, 0.0);

    for (int i = 0; i <= m_params.fftSize/2; i++) {
        mags[i] = cframe[i*2] * cframe[i*2] + cframe[i*2+1] * cframe[i*2+1];
    }

    return processMags(mags);
}

feature_t
FeatureExtractor::processMags(const vector<float> &mags)
{
    feature_t frame(m_featureSize, 0.0);

    if (!m_params.useChromaFrequencyMap &&
        (m_params.referenceFrequency != 440.)) {

        // See comment in makeStandardFrequencyMap above
        vector<float> scaled = scaleMags(mags);

        for (int i = 0; i <= m_params.fftSize/2; i++) {
            int index = m_freqMap[i];
            if (index >= 0) {
                frame[index] += scaled[i];
            }
        }

    } else {
        for (int i = 0; i <= m_params.fftSize/2; i++) {
            int index = m_freqMap[i];
            if (index >= 0) {
                frame[index] += mags[i];
            }
        }
    }

    return frame;
}

vector<float>
FeatureExtractor::scaleMags(const vector<float> &mags)
{
    // Scale the pitch content in the given magnitude spectrum to
    // accommodate a difference in tuning frequency (between the 440Hz
    // reference and the actual tuning frequency of the input audio).
    // We only do this when not using chroma features -- see the
    // comment in makeStandardFrequencyMap() above.

    if (m_params.useChromaFrequencyMap) return mags;

    double ratio = 440.f / m_params.referenceFrequency;

    int n = static_cast<int>(mags.size());

    vector<float> scaled(n, 0.0);

    for (int target = 0; target < n; ++target) {

        double source = target / ratio;

        int lower = int(source);
        int higher = lower + 1;

        double lowerProp = higher - source;
        double higherProp = source - lower;

        double value = 0.0;
        if (lower >= 0 && lower < n) {
            value += lowerProp * mags[lower];
        }
        if (higher >= 0 && higher < n) {
            value += higherProp * mags[higher];
        }

        scaled[target] = float(value);
    }

    return scaled;
}

