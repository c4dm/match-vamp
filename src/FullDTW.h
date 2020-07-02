/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright (c) 2007-2020 Simon Dixon, Chris Cannam, and Queen Mary
    University of London, Copyright (c) 2014-2015 Tido GmbH.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef FULL_DTW_H
#define FULL_DTW_H

#include <vector>
#include <functional>

#include "DistanceMetric.h"
#include "MatchTypes.h"

/** Represents a full-matrix dynamic time warping aligner. Unlike the
 *  combination of Matcher and Finder, which implement a "live" online
 *  DTW that stores only a part of the matrix, FullDTW works offline,
 *  requiring all features up-front and storing the whole cost
 *  matrix. This is far simpler than Matcher/Finder, can calculate
 *  subsequence matches, and always returns the optimal path, but it
 *  is also much more expensive - O(mn) rather than O(m+n) in both
 *  time and space.
 */
class FullDTW
{
public:
    struct Parameters {

        Parameters(double hopTime_) :
            hopTime(hopTime_),
            diagonalWeight(1.0),
            subsequence(false)
        {}

        /** Spacing of audio frames (determines the amount of overlap
         *  or skip between frames). This value is expressed in
         *  seconds.
         */
        double hopTime;

        /** Weight applied to cost of diagonal step relative to
         *  horizontal or vertical step. The default of 1.0 is normal
         *  for cases that are not expected to differ wildly in speed.
         */
        double diagonalWeight;

        /** Whether this is a subsequence matcher. The FullDTW aligns
         *  two sequences, and by default it assumes that they are
         *  anchored to one another at both ends (i.e. they both span
         *  the same material). If subsequence is true, then it
         *  instead assumes that the second sequence is intended to
         *  match some subsequence of the first, and it returns a
         *  match against the best-matching subsequence rather than
         *  the whole of the first sequence.
         */
        bool subsequence;
    };
    
    /** Constructor for FullDTW.
     * 
     *  A FullDTW expects to be provided with two sequences of feature
     *  vectors calculated by some external code (for example, a
     *  FeatureExtractor) which it will align.
     */
    FullDTW(Parameters params, DistanceMetric::Parameters dparams);

    /**
     * Align the sequence s2 against the sequence s1, returning the
     * index into s1 for each element in s2. If the subsequence
     * parameter was set on construction, then the alignment is
     * against the best-matching subsequence of s1; otherwise it is
     * against the whole of s1.
     */
    std::vector<size_t> align(const featureseq_t &s1,
                              const featureseq_t &s2);

private:
    Parameters m_params;
    DistanceMetric m_metric;

    struct CostOption {
        bool present;
        pathcost_t cost;
    };

    pathcost_t choose(CostOption x, CostOption y, CostOption d);
    pathcostmat_t costSequences(const featureseq_t &s1, const featureseq_t &s2);
};

#endif
