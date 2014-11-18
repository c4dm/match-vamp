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

#include "MatchFeatureFeeder.h"

using std::vector;

MatchFeatureFeeder::MatchFeatureFeeder(Matcher *m1, Matcher *m2) :
    pm1(m1), pm2(m2)
{
    finder = new Finder(m1, m2);
}

MatchFeatureFeeder::~MatchFeatureFeeder()
{
    delete finder;
}

void
MatchFeatureFeeder::feed(vector<double> f1, vector<double> f2)
{
    // We maintain two FIFO queues of feature vectors, one per input
    // stream.  When the match-feeder function is entered, it knows
    // that it has at least one feature in each queue.  It loops,
    // processing up to one feature per matcher, until a queue is
    // empty.  Then it returns, to be called again with more data.

    if (!f1.empty()) {
        q1.push(f1);
    }
    
    if (!f2.empty()) {
        q2.push(f2);
    }

    while (!q1.empty() && !q2.empty()) {
        feedBlock();
    }
}

void
MatchFeatureFeeder::finish()
{
    while (!q1.empty() || !q2.empty()) {
        feedBlock();
    }
}

void
MatchFeatureFeeder::feedBlock()
{
    if (q1.empty()) { // ended
        feed2();
    } else if (q2.empty()) { // ended
        feed1();
    } else if (pm1->m_frameCount < pm1->m_blockSize) {		// fill initial block
        feed1();
        feed2();
    } else if (pm1->m_runCount >= pm1->m_params.maxRunCount) {  // slope constraints
        feed2();
    } else if (pm2->m_runCount >= pm2->m_params.maxRunCount) {
        feed1();
    } else {
        switch (finder->getExpandDirection
                (pm1->m_frameCount-1, pm2->m_frameCount-1)) {
        case ADVANCE_THIS:
            feed1();
            break;
        case ADVANCE_OTHER:
            feed2();
            break;
        case ADVANCE_BOTH:
            feed1();
            feed2();
            break;
        }
    }
}

void
MatchFeatureFeeder::feed1()
{
    pm1->consumeFeatureVector(q1.front());
    q1.pop();
}

void
MatchFeatureFeeder::feed2()
{
    pm2->consumeFeatureVector(q2.front());
    q2.pop();
}

