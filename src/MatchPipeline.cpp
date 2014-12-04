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

#include "MatchPipeline.h"

MatchPipeline::MatchPipeline(FeatureExtractor::Parameters feParams,
			     FeatureConditioner::Parameters fcParams,
			     Matcher::Parameters matchParams) :
    m_fe1(feParams),
    m_fe2(feParams),
    m_fc1(fcParams),
    m_fc2(fcParams),
    m_pm1(matchParams, 0),
    m_pm2(matchParams, &m_pm1),
    m_feeder(&m_pm1, &m_pm2),
    m_lastFrameIn1(0),
    m_lastFrameIn2(0),
    m_frameNo(0)
{
    m_pm1.setOtherMatcher(&m_pm2);
}

MatchPipeline::~MatchPipeline()
{
}

void
MatchPipeline::feedFrequencyDomainAudio(const float *arr1, const float *arr2)
{
    feedFeatures(m_fe1.process(arr1), m_fe2.process(arr2));
}

void
MatchPipeline::feedFeatures(const vector<double> &f1, const vector<double> &f2)
{
    feedConditionedFeatures(m_fc1.process(f1), m_fc2.process(f2));
}

void
MatchPipeline::feedConditionedFeatures(const vector<double> &f1, const vector<double> &f2)
{
    m_feeder.feed(f1, f2);

    if (aboveThreshold(f1)) m_lastFrameIn1 = m_frameNo;
    if (aboveThreshold(f2)) m_lastFrameIn2 = m_frameNo;

    ++m_frameNo;
}

bool
MatchPipeline::aboveThreshold(const vector<double> &f)
{
    double threshold = 1e-4f;
    double sum = 0.f;
    for (int i = 0; i < int(f.size()); ++i) {
        sum += f[i] * f[i];
    }
    return (sum >= threshold);
}

void
MatchPipeline::finish()
{
    m_feeder.finish();
    getFinder()->setDurations(m_lastFrameIn1, m_lastFrameIn2);
}

Finder *
MatchPipeline::getFinder()
{
    return m_feeder.getFinder();
}




