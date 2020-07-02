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

#include "MatchFeatureFeeder.h"

using std::vector;
using std::cerr;
using std::endl;

MatchFeatureFeeder::MatchFeatureFeeder(Matcher *m1, Matcher *m2) :
    m_pm1(m1),
    m_pm2(m2),
    m_finder(m_pm1)
{
}

MatchFeatureFeeder::~MatchFeatureFeeder()
{
}

void
MatchFeatureFeeder::feed(feature_t f1, feature_t f2)
{
    // We maintain two FIFO queues of feature vectors, one per input
    // stream.  When the match-feeder function is entered, it knows
    // that it has at least one feature in each queue.  It loops,
    // processing up to one feature per matcher, until a queue is
    // empty.  Then it returns, to be called again with more data.

    if (!f1.empty()) {
        m_q1.push(f1);
    }
    
    if (!f2.empty()) {
        m_q2.push(f2);
    }

    while (!m_q1.empty() && !m_q2.empty()) {
        feedBlock();
    }
}

int
MatchFeatureFeeder::getEstimatedReferenceFrame()
{
    if (m_pm1->getFrameCount() == 0 || m_pm2->getFrameCount() == 0) {
        return 0;
    }
    int bestRow = 0;
    normpathcost_t bestCost = 0;
    if (!m_finder.getBestColCost(m_pm2->getFrameCount()-1, bestRow, bestCost)) {
        return -1;
    } else {
        return bestRow;
    }
}

void
MatchFeatureFeeder::finish()
{
    while (!m_q1.empty() || !m_q2.empty()) {
        feedBlock();
    }

//    cerr << "MatchFeatureFeeder::finish: have " << m_pm1->getFrameCount()
//         << " reference and " << m_pm2->getFrameCount() << " other frames"
//         << endl;
}

void
MatchFeatureFeeder::feedBlock()
{
    if (m_q1.empty()) { // ended
        feed2();
    } else if (m_q2.empty()) { // ended
        feed1();
    } else if (m_pm1->isFillingInitialBlock()) {
        feed1();
        feed2();
    } else if (m_pm1->isOverrunning()) { // slope constraints
        feed2();
    } else if (m_pm2->isOverrunning()) {
        feed1();
    } else {
        switch (m_finder.getExpandDirection()) {
        case AdvanceThis:
            feed1();
            break;
        case AdvanceOther:
            feed2();
            break;
        case AdvanceBoth:
            feed1();
            feed2();
            break;
        case AdvanceNone:
            cerr << "m_finder says AdvanceNone!" << endl;
            break;
        }
    }

    m_fpx.push_back(m_pm2->getFrameCount());
    m_fpy.push_back(m_pm1->getFrameCount());
}

void
MatchFeatureFeeder::feed1()
{
    m_pm1->consumeFeatureVector(m_q1.front());
    m_q1.pop();
}

void
MatchFeatureFeeder::feed2()
{
    m_pm2->consumeFeatureVector(m_q2.front());
    m_q2.pop();
}

