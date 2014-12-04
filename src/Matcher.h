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

#include "DistanceMetric.h"

using std::vector;
using std::string;
using std::cerr;
using std::endl;

/** Represents an audio feature stream that can be matched to another
 *  audio stream of the same piece of music.  The matching algorithm
 *  uses dynamic time warping.
 */
class Matcher
{
public:
    enum Advance {
        AdvanceNone,
        AdvanceBoth,
        AdvanceThis,
        AdvanceOther
    };
    static string advanceToString(Advance a) {
        switch (a) {
        case AdvanceNone: return "AdvanceNone";
        case AdvanceBoth: return "AdvanceBoth";
        case AdvanceThis: return "AdvanceThis";
        case AdvanceOther: return "AdvanceOther";
        }
        return "(unknown)";
    }

    struct Parameters {

        Parameters(float rate_, double hopTime_, int fftSize_) :
            sampleRate(rate_),
            distanceNorm(DistanceMetric::NormaliseDistanceToLogSum),
            hopTime(hopTime_),
            fftSize(fftSize_),
            blockTime(10.0),
            maxRunCount(3),
            diagonalWeight(2.0)
        {}

        /** Sample rate of audio */
        float sampleRate;

        /** Type of distance metric normalisation */
        DistanceMetric::DistanceNormalisation distanceNorm;

        /** Spacing of audio frames (determines the amount of overlap or
         *  skip between frames). This value is expressed in
         *  seconds.
         */
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

        /** Maximum number of frames sequentially processed by this
         *  matcher, without a frame of the other matcher being
         *  processed.
         */
        int maxRunCount;

        /** Weight applied to cost of diagonal step relative to
         *  horizontal or vertical step. The default of 2.0 means that
         *  a diagonal is not favoured over horizontal+vertical
         *  combined, which is good when maintaining gross tracking of
         *  performances that may have wildly differing speeds but
         *  which also leads to quite jaggy paths. A more typical
         *  normal DTW approach for performances with similar speeds
         *  might use 1.0 or something close to it.
         */
        float diagonalWeight;
    };

    /** Constructor for Matcher.
     * 
     *  A Matcher expects to be provided with feature vectors
     *  calculated by some external code (for example, a
     *  FeatureExtractor). Call consumeFeatureVector to provide each
     *  feature frame.
     *
     *  @param p The Matcher representing the performance with which
     *  this one is going to be matched.  Some information is shared
     *  between the two matchers (currently one possesses the distance
     *  matrix and optimal path matrix).
     */
    Matcher(Parameters parameters, Matcher *p);

    /** Destructor for Matcher.
     */
    ~Matcher();

    /** Adds a link to the Matcher object representing the performance
     *  which is going to be matched to this one.
     *
     *  @param p the Matcher representing the other performance
     */
    void setOtherMatcher(Matcher *p) {
        m_otherMatcher = p;
    }

    int getBlockSize() {
        return m_blockSize;
    }

    bool isOverrunning() {
        return m_runCount >= m_params.maxRunCount;
    }
    
    int getFrameCount() { 
        return m_frameCount;
    }

    int getOtherFrameCount() {
        return m_otherMatcher->getFrameCount();
    }

    float getDiagonalWeight() {
        return m_params.diagonalWeight;
    }
    
    /** Processes a feature vector frame, presumably calculated from
     *  audio data by some external code such as a FeatureExtractor.
     *  Calculates the distance to all frames stored in the
     *  otherMatcher and stores in the distance matrix, before
     *  updating the optimal path matrix using the dynamic time
     *  warping algorithm.
     *
     *  The supplied features must always be of the same size (within
     *  any pair of Matcher objects).
     */
    void consumeFeatureVector(std::vector<double> feature);
    
    /** Tests whether a location is in range in the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return true if the location is in range
     */
    bool isInRange(int i, int j);
    
    /** Tests whether a location is available in the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return true if the location is in range and contains a valid cost
     */
    bool isAvailable(int i, int j);

    /** Returns the valid range of frames in the other Matcher for the
     *  given frame in this Matcher's minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @return the first, last pair of frame numbers for the other
     *  Matcher. Note that the last frame is exclusive (last valid
     *  frame + 1).
     */
    std::pair<int, int> getColRange(int i);

    /** Returns the valid range of frames in this Matcher for the
     *  given frame in the other Matcher's minimum cost matrix.
     *
     *  @param i the frame number of the other Matcher
     *  @return the first, last pair of frame numbers for this
     *  Matcher. Note that the last frame is exclusive (last valid
     *  frame + 1).
     */
    std::pair<int, int> getRowRange(int i);
    
    /** Retrieves a value from the distance matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the distance metric at this location
     */
    float getDistance(int i, int j);

    /** Sets a value to the distance matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @param value the distance metric to set for this location
     */
    void setDistance(int i, int j, float distance);
    
    /** Retrieves a value from the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the cost of the minimum cost path to this location
     */
    double getPathCost(int i, int j);

    /** Sets a value and an advance direction to the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @param dir the direction from which this position is reached with
     *  minimum cost
     *  @param value the cost of the minimum cost path to set for this location
     */
    void setPathCost(int i, int j, Advance dir, double value);

    /** Retrieves an advance direction from the matrix.
     * 
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the direction from which this position is reached with
     *  minimum cost
     */
    Advance getAdvance(int i, int j);
    
protected:
    /** Create internal structures and reset. */
    void init();

    /** The distXSize value has changed: resize internal buffers. */
    void size();

    /** Updates an entry in the distance matrix and the optimal path matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @param dir the direction from which this position is reached with
     *  minimum cost
     *  @param value the cost of the minimum path except the current step
     *  @param dMN the distance cost between the two frames
     */
    void updateValue(int i, int j, Advance dir, double value, float dMN);

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

    /** A block of previously seen feature frames is stored in this
     *  structure for calculation of the distance matrix as the new
     *  frames are received.  One can think of the structure of the
     *  array as a circular buffer of vectors. */
    vector<vector<double> > m_frames;

    /** The best path cost matrix. */
    vector<vector<double> > m_bestPathCost;

    /** The distance matrix. */
    vector<vector<float> > m_distance;

    /** The advance direction matrix. */
    vector<vector<Advance> > m_advance;

    /** The bounds of each row of data in the distance, path cost, and
     * advance direction matrices.*/
    vector<int> m_first;
    vector<int> m_last;

    /** Width of distance, path cost, and advance direction matrices
     * and first and last vectors */
    int m_distXSize;

    bool m_initialised;

    DistanceMetric m_metric;
};

#endif
