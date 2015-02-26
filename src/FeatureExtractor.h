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

#include "MatchTypes.h"

/**
 * Convert frequency-domain audio frames into features suitable for
 * MATCH alignment calculation.
 *
 * The default feature is a warping of the frequency data to map FFT
 * frequency bins into feature bins. The mapping is linear (1-1) until
 * the resolution reaches 2 points per semitone, then logarithmic with
 * a semitone resolution.  e.g. for 44.1kHz sampling rate and fftSize
 * of 2048 (46ms), bin spacing is 21.5Hz, which is mapped linearly for
 * bins 0-34 (0 to 732Hz), and logarithmically for the remaining bins
 * (midi notes 79 to 127, bins 35 to 83), where all energy above note
 * 127 is mapped into the final bin.
 *
 * Alternatively a chroma mapping is also available. This produces a
 * 13-bin feature by mapping all FFT bins into bin 0 until the
 * resolution reaches 1 point per semitone, then mapping each
 * subsequent bin into its corresponding semitone in the remaining 12
 * bins (where bin 1 is C).  e.g. e.g. for 44.1kHz sampling rate and
 * fftSize of 2048 (46ms), frequencies up to 361 Hz go to bin 0,
 * subsequent frequencies to the chroma bins.
 */
class FeatureExtractor
{
public:
    struct Parameters {

        Parameters(float rate_, int fftSize_) :
            sampleRate(rate_),
            useChromaFrequencyMap(false),
            fftSize(fftSize_),
            referenceFrequency(440.0),
            minFrequency(0.),
            maxFrequency(rate_/2.)
        {}

        /** Sample rate of audio */
        float sampleRate;

        /** Flag indicating whether to use a chroma frequency map (12
         *  bins) instead of the default warped spectrogram */
        bool useChromaFrequencyMap;

        /** Size of an FFT frame in samples. Note that the data passed
         *  in is already in the frequency domain, so this expresses
         *  the size of the frame that the caller will be providing. */
        int fftSize;

        /** Frequency of concert A */
        double referenceFrequency;

        /** Minimum frequency cutoff to include in feature */
        double minFrequency;

        /** Maximum frequency cutoff to include in feature */
        double maxFrequency;
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
     * Return the feature vector size that would be returned from
     * process() with these parameters.
     */
    static int getFeatureSizeFor(Parameters params);
    
    /**
     * Process one frequency-domain audio frame (provided as real &
     * imaginary components from the FFT output). Return a feature
     * vector of size given by getFeatureSize(). Input vectors must
     * have at least params.fftSize/2+1 elements each.
     *
     * Operates by mapping the frequency bins into a part-linear
     * part-logarithmic array, unless useChromaFrequencyMap is true in
     * which case they are mapped into chroma bins.
     */
    feature_t process(const std::vector<double> &real,
                      const std::vector<double> &imag);
    
    /**
     * Process one frequency-domain audio frame (provided as real &
     * imaginary components from the FFT output). Return a feature
     * vector of size given by getFeatureSize(). Input vectors must
     * have at least params.fftSize/2+1 elements each.
     *
     * Operates by mapping the frequency bins into a part-linear
     * part-logarithmic array, unless useChromaFrequencyMap is true in
     * which case they are mapped into chroma bins.
     */
    feature_t process(const std::vector<float> &real,
                      const std::vector<float> &imag);
    
    /**
     * Process one frequency-domain audio frame, provided as a single
     * array of alternating real and imaginary components. Input array
     * must have at least 2 * (params.fftSize/2 + 1) elements.
     *
     * Operates by mapping the frequency bins into a part-linear
     * part-logarithmic array, unless useChromaFrequencyMap is true in
     * which case they are mapped into chroma bins.
     */
    feature_t process(const float *carray);
    
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

    /** A mapping function for mapping FFT bins to final frequency
     *  bins.  The mapping is linear (1-1) until the resolution
     *  reaches 2 points per semitone, then logarithmic with a
     *  semitone resolution.  e.g. for 44.1kHz sampling rate and
     *  fftSize of 2048 (46ms), bin spacing is 21.5Hz, which is mapped
     *  linearly for bins 0-34 (0 to 732Hz), and logarithmically for
     *  the remaining bins (midi notes 79 to 127, bins 35 to 83),
     *  where all energy above note 127 is mapped into the final
     *  bin.
     * 
     *  If a bin's frequency is outside the minFrequency->maxFrequency
     *  range, it will be mapped to a target bin of -1 and should be
     *  discarded.
     */
    std::vector<int> m_freqMap;

    feature_t processMags(const std::vector<float> &mags);
    std::vector<float> scaleMags(const std::vector<float> &mags);
    
    /** The size of a returned feature. */
    int m_featureSize;
};

#endif
    
