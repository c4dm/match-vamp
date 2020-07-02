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

#ifndef MATCHER_H
#define MATCHER_H

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>

#include "DistanceMetric.h"
#include "MatchTypes.h"

/** Represents an audio feature stream that can be matched to another
 *  audio stream of the same piece of music.  The matching algorithm
 *  uses dynamic time warping.
 */
class Matcher
{
public:
    static std::string advanceToString(advance_t a) {
        switch (a) {
        case AdvanceNone: return "AdvanceNone";
        case AdvanceBoth: return "AdvanceBoth";
        case AdvanceThis: return "AdvanceThis";
        case AdvanceOther: return "AdvanceOther";
        }
        return "(unknown)";
    }

    struct Parameters {

        Parameters(double hopTime_) :
            hopTime(hopTime_),
            blockTime(10.0),
            maxRunCount(3),
            diagonalWeight(2.0)
        {}

        /** Spacing of audio frames (determines the amount of overlap
         *  or skip between frames). This value is expressed in
         *  seconds.
         */
        double hopTime;
        
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
        double diagonalWeight;
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
    Matcher(Parameters params, DistanceMetric::Parameters dparams, Matcher *p);

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

    bool isFillingInitialBlock() {
        return m_frameCount < m_blockSize;
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

    double getDiagonalWeight() {
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
    void consumeFeatureVector(const feature_t &feature);
    
    /** Tests whether a location is in range in the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return true if the location is in range
     */
    bool isInRange(int i, int j);
    
    /** Tests whether a location is available in the minimum cost
     *  matrix, that is, whether it is in range and contains a valid
     *  cost value.  Note this and its associated isRowAvailable,
     *  isColAvailable checks are more expensive than isInRange and
     *  are really intended for error checking. (If a row is in range,
     *  it should always be available.)
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return true if the location is in range and contains a valid cost
     */
    bool isAvailable(int i, int j);

    /** Tests whether any locations in the given row are available.
     */
    bool isRowAvailable(int i);

    /** Tests whether any locations in the given column are available.
     */
    bool isColAvailable(int i);

    /** Returns the valid range of columns for the given row, that is,
     *  the range of frames in the other Matcher for the given frame
     *  in this Matcher's minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @return the first, last pair of frame numbers for the other
     *  Matcher. Note that the last frame is exclusive (last valid
     *  frame + 1).
     */
    std::pair<int, int> getColRangeForRow(int i);

    /** Returns the valid range of rows for the given column, that is,
     *  the range of frames in this Matcher for the given frame in the
     *  other Matcher's minimum cost matrix.
     *
     *  @param i the frame number of the other Matcher
     *  @return the first, last pair of frame numbers for this
     *  Matcher. Note that the last frame is exclusive (last valid
     *  frame + 1).
     */
    std::pair<int, int> getRowRangeForCol(int i);
    
    /** Retrieves a value from the distance matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the distance metric at this location
     */
    distance_t getDistance(int i, int j);

    /** Sets a value to the distance matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @param value the distance metric to set for this location
     */
    void setDistance(int i, int j, distance_t distance);
    
    /** Retrieves a value from the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the cost of the minimum cost path to this location
     */
    pathcost_t getPathCost(int i, int j);

    /** Sets a value and an advance direction to the minimum cost matrix.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @param dir the direction from which this position is reached with
     *  minimum cost
     *  @param value the cost of the minimum cost path to set for this location
     */
    void setPathCost(int i, int j, advance_t dir, pathcost_t value);
    
    /** Retrieves a value from the minimum cost matrix, normalised for
     *  path length.
     *
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the cost of the minimum cost path to this location,
     *     normalised by the Manhattan distance from 0,0 to i,j
     */
    normpathcost_t getNormalisedPathCost(int i, int j);

    /** Retrieves an advance direction from the matrix.
     * 
     *  @param i the frame number of this Matcher
     *  @param j the frame number of the other Matcher
     *  @return the direction from which this position is reached with
     *  minimum cost
     */
    advance_t getAdvance(int i, int j);

    struct MemoryStats {
        double features_k;
        double pathcosts_k;
        double distances_k;
        double advances_k;
        double total_k() const {
            return features_k + pathcosts_k + distances_k + advances_k;
        }
    };
    
    /** Obtain some stats about memory consumption.
     */
    MemoryStats getMemoryStats() const;
    
    /** Print some stats about memory consumption etc to stderr.
     */
    void printStats();
    
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
    void updateValue(int i, int j, advance_t dir, pathcost_t value, distance_t dMN);

    void calcAdvance();

    /**
     * Add the given distance increment to the given path cost, and
     * return the result clipped (if necessary) at MaxPathCost.
     */
    pathcost_t addToCost(pathcost_t cost, pathcost_t increment);
    
    /** Points to the other performance with which this one is being
     *  compared.  (See <code>m_firstPM</code>)
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
    featureseq_t m_features;

    /** The best path cost matrix. */
    pathcostmat_t m_bestPathCost;

    /** The distance matrix. */
    distancemat_t m_distance;

    /** The advance direction matrix. */
    advancemat_t m_advance;

    /** The bounds of each row of data in the distance, path cost, and
     * advance direction matrices.*/
    std::vector<int> m_first;
    std::vector<int> m_last;
    
    /** Width of distance, path cost, and advance direction matrices
     * and first and last vectors */
    int m_distXSize;

    bool m_initialised;

    DistanceMetric m_metric;
};

inline Matcher::MemoryStats operator+(const Matcher::MemoryStats &a,
                                      const Matcher::MemoryStats &b)
{
    Matcher::MemoryStats m = a;
    m.features_k += b.features_k;
    m.pathcosts_k += b.pathcosts_k;
    m.distances_k += b.distances_k;
    m.advances_k += b.advances_k;
    return m;
}

#endif
