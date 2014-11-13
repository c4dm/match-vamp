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

#include "MatchFeeder.h"

using std::vector;

MatchFeeder::MatchFeeder(Matcher *m1, Matcher *m2) :
    pm1(m1), pm2(m2)
{
    fftSize = m1->m_params.fftSize;
    finder = new Finder(m1, m2);
    reBuffer = new double[fftSize/2+1];
    imBuffer = new double[fftSize/2+1];
}

MatchFeeder::~MatchFeeder()
{
    delete[] imBuffer;
    delete[] reBuffer;
    while (!q1.empty()) {
        delete[] q1.front();
        q1.pop();
    }
    while (!q2.empty()) {
        delete[] q2.front();
        q2.pop();
    }
    delete finder;
}

void
MatchFeeder::feed(const float *const *input)
{
    // We maintain two FIFO queues of audio data frame block pointers,
    // one per input stream.  When the match-feeder function is
    // entered, it knows that it has at least one block in each queue.
    // It loops, processing up to one block per matcher, until a queue
    // is empty.  Then it returns, to be called again with more data.

    prepare(input);

    while (!q1.empty() && !q2.empty()) {
//        std::cerr << "MatchFeeder::feed: q1 " << q1.size() << " q2 " << q2.size() << std::endl;
        (void)feedBlock();
    }
}

MatchFeeder::Features
MatchFeeder::feedAndGetFeatures(const float *const *input)
{
    prepare(input);

    Features all;

    while (!q1.empty() && !q2.empty()) {
        Features ff = feedBlock();
        all.f1.insert(all.f1.end(), ff.f1.begin(), ff.f1.end());
        all.f2.insert(all.f2.end(), ff.f2.begin(), ff.f2.end());
    }

    return all;
}

void
MatchFeeder::prepare(const float *const *input)
{
    float *block = new float[fftSize+2];
    for (size_t i = 0; i < fftSize+2; ++i) {
        block[i] = input[0][i];
    }
    q1.push(block);

    block = new float[fftSize+2];
    for (size_t i = 0; i < fftSize+2; ++i) {
        block[i] = input[1][i];
    }
    q2.push(block);
}

MatchFeeder::Features
MatchFeeder::feedBlock()
{
    Features ff;
    vector<double> f1, f2;

    if (pm1->m_frameCount < pm1->m_blockSize) {		// fill initial block
//        std::cerr << "feeding initial block" << std::endl;
        f1 = feed1();
        f2 = feed2();
    }
//!!!    } else if (pm1->atEnd) {
//        feed2();
//!!!    } else if (pm2->atEnd)
//        feed1();
    else if (pm1->m_runCount >= pm1->m_params.maxRunCount) {  // slope constraints
//        std::cerr << "pm1 too slopey" << std::endl;
        f2 = feed2();
    } else if (pm2->m_runCount >= pm2->m_params.maxRunCount) {
//        std::cerr << "pm2 too slopey" << std::endl;
        f1 = feed1();
    } else {
        switch (finder->getExpandDirection
                (pm1->m_frameCount-1, pm2->m_frameCount-1)) {
        case ADVANCE_THIS:
//            std::cerr << "finder says ADVANCE_THIS" << std::endl;
            f1 = feed1();
            break;
        case ADVANCE_OTHER:
//            std::cerr << "finder says ADVANCE_OTHER" << std::endl;
            f2 = feed2();
            break;
        case ADVANCE_BOTH:
//            std::cerr << "finder says ADVANCE_BOTH" << std::endl;
            f1 = feed1();
            f2 = feed2();
            break;
        }
    }

    if (!f1.empty()) ff.f1.push_back(f1);
    if (!f2.empty()) ff.f2.push_back(f2);
    return ff;
}

vector<double>
MatchFeeder::feed1()
{
//    std::cerr << "feed1" << std::endl;
    float *block = q1.front();
    q1.pop();
    for (size_t i = 0; i <= fftSize/2; ++i) {
        reBuffer[i] = block[i*2];
    }
    for (size_t i = 0; i <= fftSize/2; ++i) {
        imBuffer[i] = block[i*2+1];
    }
    delete[] block;
    return pm1->consumeFrame(reBuffer, imBuffer);
}

vector<double>
MatchFeeder::feed2()
{
//    std::cerr << "feed2" << std::endl;
    float *block = q2.front();
    q2.pop();
    for (size_t i = 0; i <= fftSize/2; ++i) {
        reBuffer[i] = block[i*2];
    }
    for (size_t i = 0; i <= fftSize/2; ++i) {
        imBuffer[i] = block[i*2+1];
    }
    delete[] block;
    return pm2->consumeFrame(reBuffer, imBuffer);
}

