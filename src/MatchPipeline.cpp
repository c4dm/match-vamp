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

//#define DEBUG_MATCH_PIPELINE 1

MatchPipeline::MatchPipeline(FeatureExtractor::Parameters feParams,
			     FeatureConditioner::Parameters fcParams,
                             DistanceMetric::Parameters dParams,
			     Matcher::Parameters matchParams,
                             double secondReferenceFrequency) :
    m_fe1(feParams),
    m_fe2(feParams),
    m_fc1(fcParams),
    m_fc2(fcParams),
    m_pm1(matchParams, dParams, 0),
    m_pm2(matchParams, dParams, &m_pm1),
    m_feeder(&m_pm1, &m_pm2),
    m_lastFrameIn1(0),
    m_lastFrameIn2(0),
    m_frameNo(0)
{
    if (secondReferenceFrequency != 0.0) {
        feParams.referenceFrequency = secondReferenceFrequency;
        m_fe2 = FeatureExtractor(feParams);
    }
    
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
    m_f1 = f1;
    m_f2 = f2;

#ifdef DEBUG_MATCH_PIPELINE
    if (m_lastFrameIn1 == 1) {
        cerr << "features 1 -> ";
        for (int i = 0; i < (int) m_f1.size(); ++i) {
            cerr << m_f1[i] << " ";
        }
        cerr << endl;
    }
#endif
    
    feedConditionedFeatures(m_fc1.process(f1), m_fc2.process(f2));
}

void
MatchPipeline::feedConditionedFeatures(const vector<double> &c1, const vector<double> &c2)
{
    m_c1 = c1;
    m_c2 = c2;

#ifdef DEBUG_MATCH_PIPELINE
    if (m_lastFrameIn1 == 1) {
        cerr << "conditioned features 1 -> ";
        for (int i = 0; i < (int) m_c1.size(); ++i) {
            cerr << m_c1[i] << " ";
        }
        cerr << endl;
    }
#endif
    
    m_feeder.feed(c1, c2);

    if (aboveThreshold(c1)) m_lastFrameIn1 = m_frameNo;
    if (aboveThreshold(c2)) m_lastFrameIn2 = m_frameNo;

#ifdef DEBUG_MATCH_PIPELINE
    cerr << "last frames are " << m_lastFrameIn1 << ", " << m_lastFrameIn2
         << endl;
#endif
    
    ++m_frameNo;
}

void
MatchPipeline::extractFeatures(vector<double> &f1, vector<double> &f2)
{
    f1 = m_f1;
    f2 = m_f2;
}

void
MatchPipeline::extractConditionedFeatures(vector<double> &c1, vector<double> &c2)
{
    c1 = m_c1;
    c2 = m_c2;
}

bool
MatchPipeline::aboveThreshold(const vector<double> &f)
{
    // This threshold is used only to determine when either of the
    // input streams has ended -- the last frame for a stream is
    // considered to be the last one that was above the
    // threshold. This is different from the silence threshold in
    // FeatureConditioner.
    double threshold = 1e-4f;
    double sum = 0.f;
    for (int i = 0; i < int(f.size()); ++i) {
        sum += f[i] * f[i];
    }
#ifdef DEBUG_MATCH_PIPELINE
    cerr << "aboveThreshold: sum " << sum << ", threshold " << threshold
         << ", returning " << (sum >= threshold) << endl;
#endif
    return (sum >= threshold);
}

void
MatchPipeline::finish()
{
    m_feeder.finish();
    getFinder()->setDurations(m_lastFrameIn1, m_lastFrameIn2);
}

MatchFeatureFeeder *
MatchPipeline::getFeeder()
{
    return &m_feeder;
}

Finder *
MatchPipeline::getFinder()
{
    return m_feeder.getFinder();
}




