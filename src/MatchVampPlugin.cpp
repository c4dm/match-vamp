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

#include "MatchVampPlugin.h"

#include "Matcher.h"
#include "MatchFeatureFeeder.h"
#include "FeatureExtractor.h"
#include "Path.h"

#include <vamp/vamp.h>
#include <vamp-sdk/PluginAdapter.h>
#include <vamp-sdk/RealTime.h>

#include <vector>
#include <algorithm>

//static int extant = 0;

#ifdef _WIN32
HANDLE
MatchVampPlugin::m_serialisingMutex;
#else
pthread_mutex_t 
MatchVampPlugin::m_serialisingMutex;
#endif

bool
MatchVampPlugin::m_serialisingMutexInitialised = false;

// We want to ensure our freq map / crossover bin in Matcher.cpp are
// always valid with a fixed FFT length in seconds, so must reject low
// sample rates
static float sampleRateMin = 5000.f;

static float defaultStepTime = 0.020f;

MatchVampPlugin::MatchVampPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(int(inputSampleRate * defaultStepTime + 0.001)),
    m_stepTime(defaultStepTime),
    m_blockSize(2048),
    m_serialise(false),
    m_begin(true),
    m_locked(false),
    m_smooth(false),
    m_frameNo(0),
    m_params(defaultStepTime),
    m_defaultParams(defaultStepTime),
    m_feParams(inputSampleRate, m_blockSize),
    m_defaultFeParams(inputSampleRate, m_blockSize),
    m_secondReferenceFrequency(m_defaultFeParams.referenceFrequency),
    m_fcParams(),
    m_defaultFcParams(),
    m_dParams(),
    m_defaultDParams()
{
    if (inputSampleRate < sampleRateMin) {
        std::cerr << "MatchVampPlugin::MatchVampPlugin: input sample rate "
                  << inputSampleRate << " < min supported rate "
                  << sampleRateMin << ", plugin will refuse to initialise"
                  << std::endl;
    }

    if (!m_serialisingMutexInitialised) {
        m_serialisingMutexInitialised = true;
#ifdef _WIN32
        m_serialisingMutex = CreateMutex(NULL, FALSE, NULL);
#else
        pthread_mutex_init(&m_serialisingMutex, 0);
#endif
    }

    m_pipeline = 0;
//    std::cerr << "MatchVampPlugin::MatchVampPlugin(" << this << "): extant = " << ++extant << std::endl;
}

MatchVampPlugin::~MatchVampPlugin()
{
//    std::cerr << "MatchVampPlugin::~MatchVampPlugin(" << this << "): extant = " << --extant << std::endl;

    delete m_pipeline;

    if (m_locked) {
#ifdef _WIN32
        ReleaseMutex(m_serialisingMutex);
#else
        pthread_mutex_unlock(&m_serialisingMutex);
#endif
        m_locked = false;
    }
}

string
MatchVampPlugin::getIdentifier() const
{
    return "match";
}

string
MatchVampPlugin::getName() const
{
    return "Match Performance Aligner";
}

string
MatchVampPlugin::getDescription() const
{
    return "Calculate alignment between two performances in separate channel inputs";
}

string
MatchVampPlugin::getMaker() const
{
    return "Simon Dixon (plugin by Chris Cannam)";
}

int
MatchVampPlugin::getPluginVersion() const
{
    return 3;
}

string
MatchVampPlugin::getCopyright() const
{
    return "GPL";
}

MatchVampPlugin::ParameterList
MatchVampPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor desc;

    desc.identifier = "freq1";
    desc.name = "Tuning frequency of first input";
    desc.description = "Tuning frequency (concert A) for the reference audio.";
    desc.minValue = 220.0;
    desc.maxValue = 880.0;
    desc.defaultValue = float(m_defaultFeParams.referenceFrequency);
    desc.isQuantized = false;
    desc.unit = "Hz";
    list.push_back(desc);

    desc.identifier = "freq2";
    desc.name = "Tuning frequency of second input";
    desc.description = "Tuning frequency (concert A) for the other audio.";
    desc.minValue = 220.0;
    desc.maxValue = 880.0;
    desc.defaultValue = float(m_defaultFeParams.referenceFrequency);
    desc.isQuantized = false;
    desc.unit = "Hz";
    list.push_back(desc);

    desc.identifier = "minfreq";
    desc.name = "Minimum frequency";
    desc.description = "Minimum frequency to include in features.";
    desc.minValue = 0.0;
    desc.maxValue = float(m_inputSampleRate / 4.f);
    desc.defaultValue = float(m_defaultFeParams.minFrequency);
    desc.isQuantized = false;
    desc.unit = "Hz";
    list.push_back(desc);

    desc.identifier = "maxfreq";
    desc.name = "Maximum frequency";
    desc.description = "Maximum frequency to include in features.";
    desc.minValue = 1000.0;
    desc.maxValue = float(m_inputSampleRate / 2.f);
    desc.defaultValue = float(m_defaultFeParams.maxFrequency);
    desc.isQuantized = false;
    desc.unit = "Hz";
    list.push_back(desc);
    
    desc.unit = "";
    
    desc.identifier = "usechroma";
    desc.name = "Feature type";
    desc.description = "Whether to use warped spectrogram or chroma frequency map";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = m_defaultFeParams.useChromaFrequencyMap ? 1 : 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    desc.valueNames.push_back("Spectral");
    desc.valueNames.push_back("Chroma");
    list.push_back(desc);

    desc.valueNames.clear();

    desc.identifier = "usespecdiff";
    desc.name = "Use feature difference";
    desc.description = "Whether to use half-wave rectified feature-to-feature difference instead of straight spectral or chroma feature";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = float(m_defaultFcParams.order);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "framenorm";
    desc.name = "Frame normalisation";
    desc.description = "Type of normalisation to use for features";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = float(m_defaultFcParams.norm);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    desc.valueNames.push_back("None");
    desc.valueNames.push_back("Sum to 1");
    desc.valueNames.push_back("Long-term average");
    list.push_back(desc);
    desc.valueNames.clear();
    desc.defaultValue = float(m_defaultFcParams.silenceThreshold);

    desc.identifier = "metric";
    desc.name = "Distance metric";
    desc.description = "Metric for distance calculations.";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = float(m_defaultDParams.metric);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    desc.valueNames.push_back("Manhattan");
    desc.valueNames.push_back("Euclidean");
    desc.valueNames.push_back("Cosine");
    list.push_back(desc);
    desc.valueNames.clear();

    desc.identifier = "distnorm";
    desc.name = "Distance normalisation";
    desc.description = "Type of normalisation to use for distance metric";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = float(m_defaultDParams.norm);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    desc.valueNames.push_back("None");
    desc.valueNames.push_back("Sum of frames");
    desc.valueNames.push_back("Log sum of frames");
    list.push_back(desc);
    desc.valueNames.clear();

    desc.identifier = "silencethreshold";
    desc.name = "Silence threshold";
    desc.description = "Total frame energy threshold below which a feature will be regarded as silent";
    desc.minValue = 0;
    desc.maxValue = 0.1f;
    desc.defaultValue = float(m_defaultFcParams.silenceThreshold);
    desc.isQuantized = false;
    list.push_back(desc);

    desc.identifier = "noise";
    desc.name = "Add noise";
    desc.description = "Whether to mix in a small constant white noise term when calculating feature distance. This can improve alignment against sources containing cleanly synthesised audio.";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = float(m_defaultDParams.noise);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);
    
    desc.identifier = "gradientlimit";
    desc.name = "Gradient limit";
    desc.description = "Limit of number of frames that will be accepted from one source without a frame from the other source being accepted";
    desc.minValue = 1;
    desc.maxValue = 10;
    desc.defaultValue = float(m_defaultParams.maxRunCount);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "zonewidth";
    desc.name = "Search zone width";
    desc.description = "Width of the search zone (error margin) either side of the ongoing match position, in seconds";
    desc.minValue = 1;
    desc.maxValue = 60;
    desc.defaultValue = float(m_defaultParams.blockTime);
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.unit = "s";
    list.push_back(desc);

    desc.identifier = "diagonalweight";
    desc.name = "Diagonal weight";
    desc.description = "Weight applied to cost of diagonal step relative to horizontal or vertical step. The default of 2.0 is good for gross tracking of quite different performances; closer to 1.0 produces a smoother path for performances more similar in tempo";
    desc.minValue = 1.0;
    desc.maxValue = 2.0;
    desc.defaultValue = float(m_defaultParams.diagonalWeight);
    desc.isQuantized = false;
    desc.unit = "";
    list.push_back(desc);
    
    desc.identifier = "smooth";
    desc.name = "Use path smoothing";
    desc.description = "Smooth the path by replacing steps with diagonals. (This was enabled by default in earlier versions of the MATCH plugin, but the default now is to produce an un-smoothed path.)";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.unit = "";
    list.push_back(desc);

    desc.identifier = "serialise";
    desc.name = "Serialise plugin invocations";
    desc.description = "Reduce potential memory load at the expense of multiprocessor performance by serialising multi-threaded plugin runs";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);
    
    return list;
}

float
MatchVampPlugin::getParameter(std::string name) const
{
    if (name == "serialise") {
        return m_serialise ? 1.0 : 0.0; 
    } else if (name == "framenorm") {
        return float(m_fcParams.norm);
    } else if (name == "distnorm") {
        return float(m_dParams.norm);
    } else if (name == "usespecdiff") {
        return float(m_fcParams.order);
    } else if (name == "usechroma") {
        return m_feParams.useChromaFrequencyMap ? 1.0 : 0.0;
    } else if (name == "gradientlimit") {
        return float(m_params.maxRunCount);
    } else if (name == "diagonalweight") {
        return float(m_params.diagonalWeight);
    } else if (name == "zonewidth") {
        return float(m_params.blockTime);
    } else if (name == "smooth") {
        return m_smooth ? 1.0 : 0.0;
    } else if (name == "silencethreshold") {
        return float(m_fcParams.silenceThreshold);
    } else if (name == "metric") {
        return float(m_dParams.metric);
    } else if (name == "noise") {
        return m_dParams.noise;
    } else if (name == "freq1") {
        return float(m_feParams.referenceFrequency);
    } else if (name == "freq2") {
        return float(m_secondReferenceFrequency);
    } else if (name == "minfreq") {
        return float(m_feParams.minFrequency);
    } else if (name == "maxfreq") {
        return float(m_feParams.maxFrequency);
    }
    
    return 0.0;
}

void
MatchVampPlugin::setParameter(std::string name, float value)
{
    if (name == "serialise") {
        m_serialise = (value > 0.5);
    } else if (name == "framenorm") {
        m_fcParams.norm = FeatureConditioner::Normalisation(int(value + 0.1));
    } else if (name == "distnorm") {
        m_dParams.norm = DistanceMetric::DistanceNormalisation(int(value + 0.1));
    } else if (name == "usespecdiff") {
        m_fcParams.order = FeatureConditioner::OutputOrder(int(value + 0.1));
    } else if (name == "usechroma") {
        m_feParams.useChromaFrequencyMap = (value > 0.5);
    } else if (name == "gradientlimit") {
        m_params.maxRunCount = int(value + 0.1);
    } else if (name == "diagonalweight") {
        m_params.diagonalWeight = value;
    } else if (name == "zonewidth") {
        m_params.blockTime = value;
    } else if (name == "smooth") {
        m_smooth = (value > 0.5);
    } else if (name == "silencethreshold") {
        m_fcParams.silenceThreshold = value;
    } else if (name == "metric") {
        m_dParams.metric = DistanceMetric::Metric(int(value + 0.1));
    } else if (name == "noise") {
        m_dParams.noise = DistanceMetric::NoiseAddition(int(value + 0.1));
    } else if (name == "freq1") {
        m_feParams.referenceFrequency = value;
    } else if (name == "freq2") {
        m_secondReferenceFrequency = value;
    } else if (name == "minfreq") {
        m_feParams.minFrequency = value;
    } else if (name == "maxfreq") {
        m_feParams.maxFrequency = value;
    }
}

size_t
MatchVampPlugin::getPreferredStepSize() const
{
    return int(m_inputSampleRate * defaultStepTime + 0.001);
}

size_t
MatchVampPlugin::getPreferredBlockSize() const
{
    return 2048;
}

void
MatchVampPlugin::createMatchers()
{
    m_params.hopTime = m_stepTime;
    m_feParams.fftSize = m_blockSize;

    m_pipeline = new MatchPipeline(m_feParams, m_fcParams, m_dParams, m_params,
                                   m_secondReferenceFrequency);
}

bool
MatchVampPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_inputSampleRate < sampleRateMin) {
        std::cerr << "MatchVampPlugin::MatchVampPlugin: input sample rate "
                  << m_inputSampleRate << " < min supported rate "
                  << sampleRateMin << std::endl;
        return false;
    }
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;
    if (stepSize > blockSize/2 ||
        blockSize != getPreferredBlockSize()) return false;

    m_stepSize = int(stepSize);
    m_stepTime = float(stepSize) / m_inputSampleRate;
    m_blockSize = int(blockSize);

    createMatchers();
    m_begin = true;
    m_locked = false;

    return true;
}

void
MatchVampPlugin::reset()
{
    delete m_pipeline;
    m_pipeline = 0;
    m_frameNo = 0;
    createMatchers();
    m_begin = true;
    m_locked = false;
}

MatchVampPlugin::OutputList
MatchVampPlugin::getOutputDescriptors() const
{
    OutputList list;

    float outRate = 1.0f / m_stepTime;

    OutputDescriptor desc;
    desc.identifier = "path";
    desc.name = "Path";
    desc.description = "Alignment path";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    m_pathOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "a_b";
    desc.name = "A-B Timeline";
    desc.description = "Timing in performance B corresponding to moments in performance A";
    desc.unit = "sec";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    m_abOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "b_a";
    desc.name = "B-A Timeline";
    desc.description = "Timing in performance A corresponding to moments in performance B";
    desc.unit = "sec";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    m_baOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "a_b_divergence";
    desc.name = "A-B Divergence";
    desc.description = "Difference between timings in performances A and B";
    desc.unit = "sec";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    m_abDivOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "a_b_temporatio";
    desc.name = "A-B Tempo Ratio";
    desc.description = "Ratio of tempi between performances A and B";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    m_abRatioOutNo = int(list.size());
    list.push_back(desc);

    int featureSize = FeatureExtractor(m_feParams).getFeatureSize();
    
    desc.identifier = "a_features";
    desc.name = "Raw A Features";
    desc.description = "Spectral features extracted from performance A";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = featureSize;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_aFeaturesOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "b_features";
    desc.name = "Raw B Features";
    desc.description = "Spectral features extracted from performance B";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = featureSize;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_bFeaturesOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "a_cfeatures";
    desc.name = "Conditioned A Features";
    desc.description = "Spectral features extracted from performance A, after normalisation and conditioning";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = featureSize;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_caFeaturesOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "b_cfeatures";
    desc.name = "Conditioned B Features";
    desc.description = "Spectral features extracted from performance B, after norrmalisation and conditioning";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = featureSize;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_cbFeaturesOutNo = int(list.size());
    list.push_back(desc);

    desc.identifier = "overall_cost";
    desc.name = "Overall Cost";
    desc.description = "Normalised overall path cost for the cheapest path";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = 1;
    m_overallCostOutNo = int(list.size());
    list.push_back(desc);
    
    return list;
}

MatchVampPlugin::FeatureSet
MatchVampPlugin::process(const float *const *inputBuffers,
                         Vamp::RealTime timestamp)
{
    if (m_begin) {
        if (!m_locked && m_serialise) {
            m_locked = true;
#ifdef _WIN32
            WaitForSingleObject(m_serialisingMutex, INFINITE);
#else
            pthread_mutex_lock(&m_serialisingMutex);
#endif
        }
        m_startTime = timestamp;
        m_begin = false;
    }
    
//    std::cerr << timestamp.toString();

    m_pipeline->feedFrequencyDomainAudio(inputBuffers[0], inputBuffers[1]);

    FeatureSet returnFeatures;

    feature_t f1, f2;
    m_pipeline->extractFeatures(f1, f2);

    feature_t cf1, cf2;
    m_pipeline->extractConditionedFeatures(cf1, cf2);

    Feature f;
    f.hasTimestamp = false;

    f.values.clear();
    for (auto v: f1) f.values.push_back(float(v));
    returnFeatures[m_aFeaturesOutNo].push_back(f);

    f.values.clear();
    for (auto v: f2) f.values.push_back(float(v));
    returnFeatures[m_bFeaturesOutNo].push_back(f);

    f.values.clear();
    for (auto v: cf1) f.values.push_back(float(v));
    returnFeatures[m_caFeaturesOutNo].push_back(f);

    f.values.clear();
    for (auto v: cf2) f.values.push_back(float(v));
    returnFeatures[m_cbFeaturesOutNo].push_back(f);

//    std::cerr << ".";
//    std::cerr << std::endl;

    ++m_frameNo;
    
    return returnFeatures;
}

MatchVampPlugin::FeatureSet
MatchVampPlugin::getRemainingFeatures()
{
    m_pipeline->finish();

    FeatureSet returnFeatures;
    
    std::vector<int> pathx;
    std::vector<int> pathy;
    int len = m_pipeline->retrievePath(m_smooth, pathx, pathy);

    double cost = m_pipeline->getOverallCost();
    Feature costFeature;
    costFeature.hasTimestamp = false;
    costFeature.values.push_back(float(cost));
    returnFeatures[m_overallCostOutNo].push_back(costFeature);
    
    int prevx = 0;
    int prevy = 0;

    for (int i = 0; i < len; ++i) {

        int x = pathx[i];
        int y = pathy[i];

        Vamp::RealTime xt = Vamp::RealTime::frame2RealTime
            (x * m_stepSize, int(m_inputSampleRate + 0.5));
        Vamp::RealTime yt = Vamp::RealTime::frame2RealTime
            (y * m_stepSize, int(m_inputSampleRate + 0.5));

        Feature feature;
        feature.hasTimestamp = true;
        feature.timestamp = m_startTime + xt;
        feature.values.clear();
        feature.values.push_back(float(yt.sec + double(yt.nsec)/1.0e9));
        returnFeatures[m_pathOutNo].push_back(feature);
        
        if (x != prevx) {

            feature.hasTimestamp = true;
            feature.timestamp = m_startTime + xt;
            feature.values.clear();
            feature.values.push_back(float(yt.sec + yt.msec()/1000.0));
            returnFeatures[m_abOutNo].push_back(feature);

            Vamp::RealTime diff = yt - xt;
            feature.values.clear();
            feature.values.push_back(float(diff.sec + diff.msec()/1000.0));
            returnFeatures[m_abDivOutNo].push_back(feature);

            if (i > 0) {
                int lookback = 100; //!!! arbitrary
                if (lookback > i) lookback = i;
                int xdiff = x - pathx[i-lookback];
                int ydiff = y - pathy[i-lookback];
                if (xdiff != 0 && ydiff != 0) {
                    float ratio = float(ydiff)/float(xdiff);
                    if (ratio < 8 && ratio > (1.0/8)) { //!!! just for now, since we aren't dealing properly with silence yet
                        feature.values.clear();
                        feature.values.push_back(ratio);
                        returnFeatures[m_abRatioOutNo].push_back(feature);
                    }
                }
            }
        }

        if (y != prevy) {
            feature.hasTimestamp = true;
            feature.timestamp = m_startTime + yt;
            feature.values.clear();
            feature.values.push_back(float(xt.sec + xt.msec()/1000.0));
            returnFeatures[m_baOutNo].push_back(feature);
        }

        prevx = x;
        prevy = y;
    }

    delete m_pipeline;
    m_pipeline = 0;

    if (m_locked) {
#ifdef _WIN32
        ReleaseMutex(m_serialisingMutex);
#else
        pthread_mutex_unlock(&m_serialisingMutex);
#endif
        m_locked = false;
    }

    return returnFeatures;
    

/*
    for (int i = 0; i < len; ++i) {
        std::cerr << i << ": [" << pathx[i] << "," << pathy[i] << "]" << std::endl;
    }

    std::cerr << std::endl;
    std::cerr << "File: A" << std::endl;
    std::cerr << "Marks: -1" << std::endl;
    std::cerr << "FixedPoints: true 0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "0" << std::endl;
    std::cerr << "File: B" << std::endl;
    std::cerr << "Marks: 0" << std::endl;
    std::cerr << "FixedPoints: true 0" << std::endl;
    std::cerr << "0.02" << std::endl;
    std::cerr << "0.02" << std::endl;

    std::cerr << len << std::endl;
    for (int i = 0; i < len; ++i) {
        std::cerr << pathx[i] << std::endl;
    }

    std::cerr << len << std::endl;
    for (int i = 0; i < len; ++i) {
        std::cerr << pathy[i] << std::endl;
    }
*/
}

static Vamp::PluginAdapter<MatchVampPlugin> mvpAdapter;

const VampPluginDescriptor *vampGetPluginDescriptor(unsigned int version,
                                                    unsigned int index)
{
    if (version < 1) return 0;

    switch (index) {
    case  0: return mvpAdapter.getDescriptor();
    default: return 0;
    }
}
