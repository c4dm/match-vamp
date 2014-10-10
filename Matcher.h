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

#ifndef _MATCHER_H_
#define _MATCHER_H_

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>

#define ADVANCE_THIS 1
#define ADVANCE_OTHER 2
#define ADVANCE_BOTH 3
#define MASK 0xfc


using std::vector;
using std::string;
using std::cerr;
using std::endl;

/** Represents an audio stream that can be matched to another audio
 *  stream of the same piece of music.  The matching algorithm uses
 *  dynamic time warping.  The distance metric is a Euclidean metric
 *  on the FFT data with the higher frequencies mapped onto a linear
 *  scale.
 */
class Matcher
{
public:
    enum FrameNormalisation {

        /** Do not normalise audio frames */
        NoFrameNormalisation,
        
        /** Normalise each frame of audio to have a sum of 1 */
        NormaliseFrameToSum1,
        
        /** Normalise each frame of audio by the long-term average
         *  of the summed energy */
        NormaliseFrameToLTAverage,
    };

    enum DistanceNormalisation {
            
        /** Do not normalise distance metrics */
        NoDistanceNormalisation,

        /** Normalise distance metric for pairs of audio frames by
         *  the sum of the two frames. */
        NormaliseDistanceToSum,

        /** Normalise distance metric for pairs of audio frames by
         *  the log of the sum of the frames. */
        NormaliseDistanceToLogSum,
    };

    struct Parameters {

        Parameters(float rate_, double hopTime_, int fftSize_) :
            sampleRate(rate_),
            frameNorm(NormaliseFrameToSum1),
            distanceNorm(NormaliseDistanceToLogSum),
            useSpectralDifference(true),
            useChromaFrequencyMap(false),
            hopTime(hopTime_),
            fftSize(fftSize_),
            blockTime(10.0),
            silenceThreshold(0.01),
            decay(0.99),
            maxRunCount(3)
        {}

        /** Sample rate of audio */
        float sampleRate;

        /** Type of audio frame normalisation */
        FrameNormalisation frameNorm;

        /** Type of distance metric normalisation */
        DistanceNormalisation distanceNorm;

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
         *  in to Matcher is already in the frequency domain, so this
         *  expresses the size of the frame that the caller will be
         *  providing.
         */
        int fftSize;

        /** The width of the search band (error margin) around the current
         *  match position, measured in seconds. Strictly speaking the
         *  width is measured backwards from the current point, since the
         *  algorithm has to work causally.
         */
        double blockTime;
        
        /** RMS level below which frame is considered silent */
        double silenceThreshold;

        /** Frame-to-frame decay factor in calculating long-term average */
        double decay;

        /** Maximum number of frames sequentially processed by this
         *  matcher, without a frame of the other matcher being
         *  processed.
         */
        int maxRunCount;
    };

protected:
    /** Points to the other performance with which this one is being
     *  compared.  The data for the distance metric and the dynamic
     *  time warping is shared between the two matchers. In the
     *  original version, only one of the two performance matchers
     *  contained the distance metric. (See <code>first</code>)
     */
    Matcher *otherMatcher;

    /** Indicates which performance is considered primary (the
     *  score). This is the performance shown on the vertical axis,
     *  and referred to as "this" in the codes for the direction of
     *  DTW steps. */
    bool firstPM;

    /** Configuration parameters */
    Parameters params;

    /** Scaling factor for distance metric; must guarantee that the
     *  final value fits in the data type used, that is, unsigned
     *  char.
     */
    double scale;

    /** Width of the search band in FFT frames (see <code>blockTime</code>) */
    int blockSize;

    /** The number of frames of audio data which have been read. */
    int frameCount;

    /** Long term average frame energy (in frequency domain
     *  representation). */
    double ltAverage;

    /** The number of frames sequentially processed by this matcher,
     *  without a frame of the other matcher being processed.
     */
    int runCount;

    /** A mapping function for mapping FFT bins to final frequency
     *  bins.  The mapping is linear (1-1) until the resolution
     *  reaches 2 points per semitone, then logarithmic with a
     *  semitone resolution.  e.g. for 44.1kHz sampling rate and
     *  fftSize of 2048 (46ms), bin spacing is 21.5Hz, which is mapped
     *  linearly for bins 0-34 (0 to 732Hz), and logarithmically for
     *  the remaining bins (midi notes 79 to 127, bins 35 to 83),
     *  where all energy above note 127 is mapped into the final
     *  bin. */
    vector<int> freqMap;

    /** The number of entries in <code>freqMap</code>. */
    int freqMapSize;

    /** The number of values in an externally-supplied feature vector,
     *  used in preference to freqMap/freqMapSize if constructed with
     *  the external feature version of the Matcher constructor. If
     *  this is zero, the internal feature extractor will be used as
     *  normal.
     */
    int externalFeatureSize;

    /** The number of values in the feature vectors actually in
     *  use. This will be externalFeatureSize if greater than zero, or
     *  freqMapSize otherwise.
     */
    int featureSize;

    /** The most recent frame; used for calculating the frame to frame
     *  spectral difference. These are therefore frequency warped but
     *  not yet normalised. */
    vector<double> prevFrame;
    vector<double> newFrame;

    /** A block of previously seen frames are stored in this structure
     *  for calculation of the distance matrix as the new frames are
     *  read in.  One can think of the structure of the array as a
     *  circular buffer of vectors.  These are the frames with all
     *  applicable processing applied (e.g. spectral difference,
     *  normalisation), unlike prevFrame and newFrame. The total
     *  energy of frames[i] is stored in totalEnergies[i]. */
    vector<vector<double> > frames;

    /** The total energy of each frame in the frames block. */ 
    vector<double> totalEnergies;

    /** The best path cost matrix. */
    int **bestPathCost;

    /** The distance matrix. */
    unsigned char **distance;

    /** The bounds of each row of data in the distance and path cost matrices.*/
    int *first;
    int *last;

    /** Height of each column in distance and bestPathCost matrices */
    int *distYSizes;

    /** Width of distance and bestPathCost matrices and first and last vectors */
    int  distXSize;

    bool initialised;

    /** Disable or enable debugging output */
    static bool silent;

public:
    /** Constructor for Matcher.
     *
     *  @param p The Matcher representing the performance with which
     *  this one is going to be matched.  Some information is shared
     *  between the two matchers (currently one possesses the distance
     *  matrix and optimal path matrix).
     */
    Matcher(Parameters parameters, Matcher *p);

    /** Constructor for Matcher using externally supplied features.
     *  A Matcher made using this constructor will not carry out its
     *  own feature extraction from frequency-domain audio data, but
     *  instead will accept arbitrary feature frames calculated by
     *  some external code.
     *
     *  @param p The Matcher representing the performance with which
     *  this one is going to be matched.  Some information is shared
     *  between the two matchers (currently one possesses the distance
     *  matrix and optimal path matrix).
     *  
     *  @param featureSize Number of values in each feature vector.
     */
    Matcher(Parameters parameters, Matcher *p, int featureSize);

    ~Matcher();

    /** For debugging, outputs information about the Matcher to
     *  standard error.
     */
    void print();

    /** Adds a link to the Matcher object representing the performance
     *  which is going to be matched to this one.
     *
     *  @param p the Matcher representing the other performance
     */
    void setOtherMatcher(Matcher *p) {
        otherMatcher = p;
    } // setOtherMatcher()

    int getFrameCount() { 
        return frameCount;
    }

    /**
     * Return the feature vector size that will be used for the given
     * parameters.
     */
    static int getFeatureSizeFor(Parameters params);

protected:
    template <typename T>
    void initVector(vector<T> &vec, int sz, T dflt = 0) {
        vec.clear();
        while ((int)vec.size() < sz) vec.push_back(dflt);
    }

    template <typename T>
    void initMatrix(vector<vector<T> > &mat, int hsz, int vsz,
                    T dflt = 0, int fillTo = -1) {
        mat.clear();
        if (fillTo < 0) fillTo = hsz;
        for (int i = 0; i < hsz; ++i) {
            mat.push_back(vector<T>());
            if (i < fillTo) {
                while ((int)mat[i].size() < vsz) {
                    mat[i].push_back(dflt);
                }
            }
        }
    }

    void init();

    void makeFreqMap();

    /** Creates a map of FFT frequency bins to comparison bins.  Where
     *  the spacing of FFT bins is less than 0.5 semitones, the
     *  mapping is one to one. Where the spacing is greater than 0.5
     *  semitones, the FFT energy is mapped into semitone-wide
     *  bins. No scaling is performed; that is the energy is summed
     *  into the comparison bins. See also consumeFrame()
     */
    void makeStandardFrequencyMap();

    void makeChromaFrequencyMap();

    /** Processes a frame of audio data by first computing the STFT
     *  with a Hamming window, then mapping the frequency bins into a
     *  part-linear part-logarithmic array, then (optionally)
     *  computing the half-wave rectified spectral difference from the
     *  previous frame, then (optionally) normalising to a sum of 1,
     *  then calculating the distance to all frames stored in the
     *  otherMatcher and storing them in the distance matrix, and
     *  finally updating the optimal path matrix using the dynamic
     *  time warping algorithm.
     *
     *  Return value is the frame (post-processed, with warping,
     *  rectification, and normalisation as appropriate).
     *
     *  The Matcher must have been constructed using the constructor
     *  without an external featureSize parameter in order to use this
     *  function. (Otherwise it will be expecting you to call
     *  consumeFeatureVector.)
     */
    std::vector<double> consumeFrame(double *reBuffer, double *imBuffer);

    /** Processes a feature vector frame (presumably calculated from
     *  audio data by some external code). As consumeFrame, except
     *  that it does not calculate a feature from audio data but
     *  instead uses the supplied feature directly.
     *
     *  The Matcher must have been constructed using the constructor
     *  that accepts an external featureSize parameter in order to
     *  use this function. The supplied feature must be of the size
     *  that was passed to the constructor.
     */
    void consumeFeatureVector(std::vector<double> feature);

    /** Calculates the Manhattan distance between two vectors, with an
     *  optional normalisation by the combined values in the
     *  vectors. Since the vectors contain energy, this could be
     *  considered as a squared Euclidean distance metric. Note that
     *  normalisation assumes the values are all non-negative.
     *
     *  @param f1 one of the vectors involved in the distance calculation
     *  @param f2 one of the vectors involved in the distance calculation
     *  @return the distance, scaled and truncated to an integer
     */
    int calcDistance(const vector<double> &f1, const vector<double> &f2);

    /** Retrieves values from the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the cost of the minimum cost path to this location
     */
    int getValue(int i, int j, bool firstAttempt);

    /** Stores entries in the distance matrix and the optimal path matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @param dir the direction from which this position is reached with
     *  minimum cost
     *  @param value the cost of the minimum path except the current step
     *  @param dMN the distance cost between the two frames
     */
    void setValue(int i, int j, int dir, int value, int dMN);

    vector<double> processFrameFromFreqData(double *, double *);
    void calcAdvance();

    friend class MatchFeeder;
    friend class MatchFeatureFeeder;
    friend class Finder;

}; // class Matcher

#endif
