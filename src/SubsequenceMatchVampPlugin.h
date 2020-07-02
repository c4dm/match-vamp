/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright (c) 2007-2020 Simon Dixon, Chris Cannam, and Queen Mary
    University of London, Copyright (c) 2014-2015 Tido GmbH.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef SUBSEQUENCE_MATCH_VAMP_PLUGIN_H
#define SUBSEQUENCE_MATCH_VAMP_PLUGIN_H

#include <vamp-sdk/Plugin.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

#include "Matcher.h"
#include "FeatureExtractor.h"
#include "FeatureConditioner.h"
#include "DistanceMetric.h"

class SubsequenceMatchVampPlugin : public Vamp::Plugin
{
public:
    SubsequenceMatchVampPlugin(float inputSampleRate);
    virtual ~SubsequenceMatchVampPlugin();

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    InputDomain getInputDomain() const { return FrequencyDomain; }

    size_t getPreferredStepSize() const;
    size_t getPreferredBlockSize() const;

    size_t getMinChannelCount() const { return 2; }
    size_t getMaxChannelCount() const { return 2; }

    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string) const;
    void setParameter(std::string, float);

    OutputList getOutputDescriptors() const;

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();

protected:
    int m_stepSize;
    float m_stepTime;
    int m_blockSize;
    int m_coarseDownsample;
    bool m_serialise;
    bool m_smooth;
    
    Matcher::Parameters m_params;
    Matcher::Parameters m_defaultParams;

    FeatureExtractor::Parameters m_feParams;
    FeatureExtractor::Parameters m_defaultFeParams;
    double m_secondReferenceFrequency;

    FeatureConditioner::Parameters m_fcParams;
    FeatureConditioner::Parameters m_defaultFcParams;

    DistanceMetric::Parameters m_dParams;
    DistanceMetric::Parameters m_defaultDParams;

    mutable int m_pathOutNo;
    mutable int m_baOutNo;
};

#endif
