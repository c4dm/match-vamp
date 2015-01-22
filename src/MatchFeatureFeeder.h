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

#ifndef _MATCH_FEATURE_FEEDER_H_
#define _MATCH_FEATURE_FEEDER_H_

#include "Matcher.h"
#include "Finder.h"

#include <queue>
#include <vector>

class MatchFeatureFeeder
{
public:
    MatchFeatureFeeder(Matcher *m1, Matcher *m2);
    ~MatchFeatureFeeder();

    /**
     * Feed the two supplied feature vectors to feeders 1 and 2
     * respectively (depending on their advance status). Matchers must
     * have been constructed using the external featureSize
     * constructor.
     *
     * f1 and f2 are normally expected to have the same number of
     * values, and that number should be the featureSize passed to the
     * constructors for both Matchers. The exception is when one input
     * ends before the other one: subsequent calls should pass a
     * feature vector as normal for the input that is still going on,
     * and an empty vector for the one that has ended.
     */
    void feed(std::vector<double> f1, std::vector<double> f2);

    /**
     * Indicate that both inputs have come to an end.
     */
    void finish();

    /**
     * Return the forward path, that is, the estimate of the
     * lowest-cost path that was generated (possibly in real-time)
     * while initially tracking the inputs. This is the path that is
     * used to determine the shape of the search zone within which the
     * eventual reverse path will be sought by the Finder.
     */
    void retrieveForwardPath(std::vector<int> &pathx, std::vector<int> &pathy) {
        pathx = m_fpx;
        pathy = m_fpy;
    }

    Finder *getFinder() { return m_finder; }

protected:
    void feedBlock();
    void feed1();
    void feed2();

    Finder *m_finder; // I own this, and it refers to m_pm1 and m_pm2
    
    Matcher *m_pm1;   // I do not own this
    Matcher *m_pm2;   // I do not own this

    std::queue<std::vector<double> > m_q1;
    std::queue<std::vector<double> > m_q2;

    vector<int> m_fpx;
    vector<int> m_fpy;
};

#endif
