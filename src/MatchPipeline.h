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

#ifndef MATCH_PIPELINE_H
#define MATCH_PIPELINE_H

#include "Matcher.h"
#include "Finder.h"
#include "FeatureExtractor.h"
#include "FeatureConditioner.h"
#include "MatchFeatureFeeder.h"

class MatchPipeline
{
public:
    /**
     * Pipeline consisting of two Matchers, two FeatureConditioners,
     * two FeatureExtractors, and a Finder. Features may be inserted
     * at any point in the pipeline.
     *
     * The pipeline goes: 
     *    Frequency-domain audio
     *      -> Features
     *          -> Conditioned features
     *              -> Matcher
     */
    MatchPipeline(FeatureExtractor::Parameters feParams,
		  FeatureConditioner::Parameters fcParams,
		  Matcher::Parameters matchParams);

    ~MatchPipeline();

    /**
     * Feed in data at the first pipeline stage. The input arrays
     * represent frames of audio from the two different sources. Each
     * is provided as a single array of alternating real and imaginary
     * components. 
     *
     * Input arrays must have at least 2 * (feParams.fftSize/2 + 1)
     * elements. The arrays will be passed to FeatureExtractor and
     * then on into the rest of the pipeline.
     */
    void feedFrequencyDomainAudio(const float *arr1, const float *arr2);

    /**
     * Feed in data at the second pipeline stage. The vectors
     * represent feature frames from two different sources. They will
     * be passed in to FeatureConditioner and then on to the rest of
     * the pipeline.
     */
    void feedFeatures(const vector<double> &f1, const vector<double> &f2);
    
    /**
     * Feed in data at the third pipeline stage. The vectors represent
     * conditioned feature frames from two different sources. They
     * will be passed to MatchFeatureFeeder for feeding to the two
     * matchers.
     */
    void feedConditionedFeatures(const vector<double> &f1, const vector<double> &f2);

    /**
     * If a frame was just fed in at the first or second pipeline
     * stage, it can be retrieved from the second stage here. That is,
     * if you provided frequency-domain audio, extractFeatures will
     * give you back the FeatureExtractor's features; if you provided
     * features, extractFeatures will simply give you those back. If
     * you provided conditioned features, the return from this
     * function is undefined.
     */
    void extractFeatures(vector<double> &f1, vector<double> &f2);

    /**
     * If a frame was just fed in at the first or second pipeline
     * stage, you can obtain the magnitudes of its features for both
     * inputs here.
     */
    void extractFeatureMagnitudes(double &mag1, double &mag2);
    
    /**
     * Retrieve the conditioned features from the third pipeline stage.
     */
    void extractConditionedFeatures(vector<double> &f1, vector<double> &f2);

    /**
     * Indicate that both inputs have come to an end.
     */
    void finish();

    MatchFeatureFeeder *getFeeder();
    Finder *getFinder();
    
private:
    FeatureExtractor m_fe1;
    FeatureExtractor m_fe2;
    FeatureConditioner m_fc1;
    FeatureConditioner m_fc2;
    Matcher m_pm1;
    Matcher m_pm2;
    MatchFeatureFeeder m_feeder;
    int m_lastFrameIn1;
    int m_lastFrameIn2;
    int m_frameNo;
    vector<double> m_f1;
    vector<double> m_f2;
    vector<double> m_c1;
    vector<double> m_c2;
    bool aboveThreshold(const vector<double> &f);
};

#endif
