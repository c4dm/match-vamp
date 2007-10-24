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

class MatchFeeder
{
public:
    MatchFeeder(Matcher *m1, Matcher *m2);
    ~MatchFeeder();

    void feed(const float *const *input);

    Finder *getFinder() { return finder; }

protected:
    void feedBlock();
    void feed1();
    void feed2();

    Finder *finder;
    Matcher *pm1;
    Matcher *pm2;

    size_t fftSize;
    double *reBuffer;
    double *imBuffer;

    std::queue<float *> q1;
    std::queue<float *> q2;
};

#endif
