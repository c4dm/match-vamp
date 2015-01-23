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
    m_pm1(m1), m_pm2(m2)
{
    m_finder = new Finder(m1);
}

MatchFeatureFeeder::~MatchFeatureFeeder()
{
    delete m_finder;
}

MatchFeatureFeeder::MatchFeatureFeeder(const MatchFeatureFeeder &other) :
    m_pm1(other.m_pm1), m_pm2(other.m_pm2)
{
    //!!! This is gross. Finder should probably not be heap allocated at all
    m_finder = new Finder(*other.m_finder);
}

MatchFeatureFeeder &
MatchFeatureFeeder::operator=(const MatchFeatureFeeder &other)
{
    m_pm1 = other.m_pm1;
    m_pm2 = other.m_pm2;
    m_finder = new Finder(*other.m_finder);
    return *this;
}

void
MatchFeatureFeeder::setMatchers(Matcher *m1, Matcher *m2)
{
    m_pm1 = m1;
    m_pm2 = m2;
    m_finder->setMatcher(m_pm1);
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
        m_q1.push(f1);
    }
    
    if (!f2.empty()) {
        m_q2.push(f2);
    }

    while (!m_q1.empty() && !m_q2.empty()) {
        feedBlock();
    }
}

void
MatchFeatureFeeder::finish()
{
    while (!m_q1.empty() || !m_q2.empty()) {
        feedBlock();
    }
}

void
MatchFeatureFeeder::feedBlock()
{
    if (m_q1.empty()) { // ended
        feed2();
    } else if (m_q2.empty()) { // ended
        feed1();
    } else if (m_pm1->getFrameCount() < m_pm1->getBlockSize()) { // fill initial block
        feed1();
        feed2();
    } else if (m_pm1->isOverrunning()) { // slope constraints
        feed2();
    } else if (m_pm2->isOverrunning()) {
        feed1();
    } else {
        switch (m_finder->getExpandDirection
                (m_pm1->getFrameCount()-1, m_pm2->getFrameCount()-1)) {
        case Matcher::AdvanceThis:
            feed1();
            break;
        case Matcher::AdvanceOther:
            feed2();
            break;
        case Matcher::AdvanceBoth:
            feed1();
            feed2();
            break;
        case Matcher::AdvanceNone:
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

