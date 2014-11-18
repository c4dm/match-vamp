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

#ifndef _MATCH_FEEDER_H_
#define _MATCH_FEEDER_H_

#include "Matcher.h"
#include "Finder.h"

#include <queue>
#include <vector>

class MatchFeeder
{
public:
    MatchFeeder(Matcher *m1, Matcher *m2);
    ~MatchFeeder();

    /**
     * Feed the two supplied channels of frequency-domain input data
     * to feeders 1 and 2 respectively, as appropriate (depending on
     * their advance status).
     */
    void feed(const float *const *input);

    /**
     * Indicate that the input has come to an end.
     */
    void finish();
    
    struct Features {
        std::vector<std::vector<double> > f1;
        std::vector<std::vector<double> > f2;
    };

    /**
     * Feed the two supplied channels of frequency-domain input data
     * to matchers 1 and 2 respectively, as appropriate (depending on
     * their advance status) and return any new feature vectors
     * calculated by the two feeders.
     */
    Features feedAndGetFeatures(const float *const *input);

    /**
     * Indicate that the input has come to an end, and return any
     * remaining features.
     */
    Features finishAndGetFeatures();
    
    Finder *getFinder() { return finder; }

protected:
    void prepare(const float *const *input);
    Features feedBlock();
    std::vector<double> feed1();
    std::vector<double> feed2();

    Finder *finder;
    Matcher *pm1;
    Matcher *pm2;

    size_t fftSize;
    double *reBuffer;
    double *imBuffer;

    std::queue<float *> q1;
    std::queue<float *> q2;

    int n;
    int lastIn1;
    int lastIn2;
};

#endif
