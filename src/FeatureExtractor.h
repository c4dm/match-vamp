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

#ifndef FEATURE_EXTRACTOR_H
#define FEATURE_EXTRACTOR_H

#include <vector>

/**
 * Convert frequency-domain audio frames into features suitable for
 * MATCH alignment calculation. The default feature is a warping of
 * the frequency data to map higher frequencies into a linear scale. A
 * chroma mapping is also available.
 *
 * Note that FeatureExtractor maintains internal frame-to-frame state:
 * use one FeatureExtractor per audio source, and construct a new one
 * for each new source.
 */
class FeatureExtractor
{
public:
    enum FrameNormalisation {

        /** Do not normalise frames */
        NoFrameNormalisation,
        
        /** Normalise each frame to have a sum of 1 */
        NormaliseFrameToSum1,
        
        /** Normalise each frame by the long-term average of the
         *  summed energy */
        NormaliseFrameToLTAverage,
    };

    struct Parameters {

        Parameters(float rate_, int fftSize_) :
            sampleRate(rate_),
            frameNorm(NormaliseFrameToSum1),
            useSpectralDifference(true),
            useChromaFrequencyMap(false),
            fftSize(fftSize_),
            silenceThreshold(0.01),
            decay(0.99)
        {}

        /** Sample rate of audio */
        float sampleRate;

        /** Type of audio frame normalisation */
        FrameNormalisation frameNorm;

        /** Flag indicating whether or not the half-wave rectified
         *  spectral difference should be used in calculating the
         *  distance metric for pairs of audio frames, instead of the
         *  straight spectrum values. */
        bool useSpectralDifference;

        /** Flag indicating whether to use a chroma frequency map (12
         *  bins) instead of the default warped spectrogram */
        bool useChromaFrequencyMap;

        /** Spacing of audio frames (determines the amount of overlap or
         *  skip between frames). This value is expressed in
         *  seconds. */
        double hopTime;

        /** Size of an FFT frame in samples. Note that the data passed
         *  in is already in the frequency domain, so this expresses
         *  the size of the frame that the caller will be providing. */
        int fftSize;

        /** RMS level below which frame is considered silent */
        double silenceThreshold;

        /** Frame-to-frame decay factor in calculating long-term average */
        double decay;
    };

    /**
     * Construct a FeatureExtractor with the given parameters.
     *
     * Note that FeatureExtractor maintains internal frame-to-frame
     * state: use one FeatureExtractor per audio source, and construct
     * a new one for each new source.
     */
    FeatureExtractor(Parameters params);

    /**
     * Return the feature vector size that will be returned from process().
     */
    int getFeatureSize() const { return m_featureSize; }
    
    /**
     * Process one frequency-domain audio frame (provided as real &
     * imaginary components from the FFT output). Return a feature
     * vector of size given by getFeatureSize().
     *
     * Operates by mapping the frequency bins into a part-linear
     * part-logarithmic array, then (optionally) computing the
     * half-wave rectified spectral difference from the previous
     * frame, then (optionally) normalising to a sum of 1.
     *
     * Return value is the frame (post-processed, with warping,
     * rectification, and normalisation as appropriate).
     */
    std::vector<double> process(const std::vector<double> &real,
                                const std::vector<double> &imag);
    
protected:
    /** Make either standard or chroma map, depending on m_params */
    void makeFreqMap();

    /** Creates a map of FFT frequency bins to comparison bins.  Where
     *  the spacing of FFT bins is less than 0.5 semitones, the
     *  mapping is one to one. Where the spacing is greater than 0.5
     *  semitones, the FFT energy is mapped into semitone-wide
     *  bins. No scaling is performed; that is the energy is summed
     *  into the comparison bins. */
    void makeStandardFrequencyMap();

    /** Creates a map of FFT frequency bins to semitone chroma bins. */
    void makeChromaFrequencyMap();

    /** Configuration parameters */
    Parameters m_params;
    
    /** Long term average frame energy (in frequency domain
     *  representation). */
    double m_ltAverage;

    /** A mapping function for mapping FFT bins to final frequency
     *  bins.  The mapping is linear (1-1) until the resolution
     *  reaches 2 points per semitone, then logarithmic with a
     *  semitone resolution.  e.g. for 44.1kHz sampling rate and
     *  fftSize of 2048 (46ms), bin spacing is 21.5Hz, which is mapped
     *  linearly for bins 0-34 (0 to 732Hz), and logarithmically for
     *  the remaining bins (midi notes 79 to 127, bins 35 to 83),
     *  where all energy above note 127 is mapped into the final
     *  bin. */
    std::vector<int> m_freqMap;

    /** The size of a returned feature. */
    int m_featureSize;

    /** The most recent frame; used for calculating the frame to frame
     *  spectral difference. This is therefore frequency warped but
     *  not yet normalised. */
    std::vector<double> m_prevFrame;
};

#endif
    
