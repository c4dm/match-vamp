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

#include "DistanceMetric.h"
#include "FeatureExtractor.h"

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
    struct Parameters {

        Parameters(float rate_, double hopTime_, int fftSize_) :
            sampleRate(rate_),
            distanceNorm(DistanceMetric::NormaliseDistanceToLogSum),
            distanceScale(90.0),
            hopTime(hopTime_),
            fftSize(fftSize_),
            blockTime(10.0),
            maxRunCount(3)
        {}

        /** Sample rate of audio */
        float sampleRate;

        /** Type of distance metric normalisation */
        DistanceMetric::DistanceNormalisation distanceNorm;

        /** Scaling factor for distance metric; must guarantee that the
         *  final value fits in the data type used, that is, unsigned
         *  char.
         */
        double distanceScale;

        /** Spacing of audio frames (determines the amount of overlap or
         *  skip between frames). This value is expressed in
         *  seconds. */
        double hopTime;
        
        /** Size of an FFT frame in samples. Note that the data passed
         *  in to Matcher is already in the frequency domain, so this
         *  expresses the size of the frame that the caller will be
         *  providing. */
        int fftSize;
        
        /** The width of the search band (error margin) around the current
         *  match position, measured in seconds. Strictly speaking the
         *  width is measured backwards from the current point, since the
         *  algorithm has to work causally.
         */
        double blockTime;

        /** Maximum number of frames sequentially processed by this
         *  matcher, without a frame of the other matcher being
         *  processed.
         */
        int maxRunCount;
    };

    /** Constructor for Matcher.
     *
     *  @param p The Matcher representing the performance with which
     *  this one is going to be matched.  Some information is shared
     *  between the two matchers (currently one possesses the distance
     *  matrix and optimal path matrix).
     */
    Matcher(Parameters parameters,
            FeatureExtractor::Parameters featureParams,
            Matcher *p);

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

    /** Adds a link to the Matcher object representing the performance
     *  which is going to be matched to this one.
     *
     *  @param p the Matcher representing the other performance
     */
    void setOtherMatcher(Matcher *p) {
        m_otherMatcher = p;
    } // setOtherMatcher()

    int getFrameCount() { 
        return m_frameCount;
    }

protected:
    /** Create internal structures and reset. */
    void init();

    /** The distXSize value has changed: resize internal buffers. */
    void size();

    /** Process a frequency-domain frame of audio data using the
     *  built-in FeatureExtractor, then calculating the distance to
     *  all frames stored in the otherMatcher and storing them in the
     *  distance matrix, and finally updating the optimal path matrix
     *  using the dynamic time warping algorithm.
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

    void calcAdvance();

    /** Points to the other performance with which this one is being
     *  compared.  The data for the distance metric and the dynamic
     *  time warping is shared between the two matchers. In the
     *  original version, only one of the two performance matchers
     *  contained the distance metric. (See <code>first</code>)
     */
    Matcher *m_otherMatcher;

    /** Indicates which performance is considered primary (the
     *  score). This is the performance shown on the vertical axis,
     *  and referred to as "this" in the codes for the direction of
     *  DTW steps. */
    bool m_firstPM;

    /** Configuration parameters */
    Parameters m_params;

    /** Width of the search band in FFT frames (see <code>blockTime</code>) */
    int m_blockSize;

    /** The number of frames of audio data which have been read. */
    int m_frameCount;

    /** The number of frames sequentially processed by this matcher,
     *  without a frame of the other matcher being processed.
     */
    int m_runCount;

    /** The number of values in a feature vector. */
    int m_featureSize;

    /** A block of previously seen frames are stored in this structure
     *  for calculation of the distance matrix as the new frames are
     *  read in.  One can think of the structure of the array as a
     *  circular buffer of vectors.  These are the frames with all
     *  applicable processing applied (e.g. spectral difference,
     *  normalisation), unlike prevFrame and newFrame. The total
     *  energy of frames[i] is stored in totalEnergies[i]. */
    vector<vector<double> > m_frames;

    /** The best path cost matrix. */
    vector<vector<int> > m_bestPathCost;

    /** The distance matrix. */
    vector<vector<unsigned char> > m_distance;

    /** The bounds of each row of data in the distance and path cost matrices.*/
    vector<int> m_first;
    vector<int> m_last;

    /** Height of each column in distance and bestPathCost matrices */
    vector<int> m_distYSizes;

    /** Width of distance and bestPathCost matrices and first and last vectors */
    int m_distXSize;

    bool m_initialised;

    FeatureExtractor m_featureExtractor;
    DistanceMetric m_metric;
    
    friend class MatchFeeder;
    friend class MatchFeatureFeeder;
    friend class Finder;

}; // class Matcher

#endif
