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
    m_smooth(true),
    m_frameNo(0),
    m_lastFrameIn1(0),
    m_lastFrameIn2(0),
    m_params(inputSampleRate, defaultStepTime, m_blockSize),
    m_defaultParams(inputSampleRate, defaultStepTime, m_blockSize),
    m_feParams(inputSampleRate, m_blockSize),
    m_defaultFeParams(inputSampleRate, m_blockSize)
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

    m_pm1 = 0;
    m_pm2 = 0;
    m_fe1 = 0;
    m_fe2 = 0;
    m_feeder = 0;
//    std::cerr << "MatchVampPlugin::MatchVampPlugin(" << this << "): extant = " << ++extant << std::endl;
}

MatchVampPlugin::~MatchVampPlugin()
{
//    std::cerr << "MatchVampPlugin::~MatchVampPlugin(" << this << "): extant = " << --extant << std::endl;

    delete m_feeder;
    delete m_fe1;
    delete m_fe2;
    delete m_pm1;
    delete m_pm2;

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
    return 2;
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

    desc.identifier = "serialise";
    desc.name = "Serialise Plugin Invocations";
    desc.description = "Reduce potential memory load at the expense of multiprocessor performance by serialising multi-threaded plugin runs";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "framenorm";
    desc.name = "Frame Normalisation";
    desc.description = "Type of normalisation to use for frequency-domain audio features";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = (int)m_defaultFeParams.frameNorm;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    desc.valueNames.push_back("None");
    desc.valueNames.push_back("Sum To 1");
    desc.valueNames.push_back("Long-Term Average");
    list.push_back(desc);
    desc.valueNames.clear();

    desc.identifier = "distnorm";
    desc.name = "Distance Normalisation";
    desc.description = "Type of normalisation to use for distance metric";
    desc.minValue = 0;
    desc.maxValue = 2;
    desc.defaultValue = (int)m_defaultParams.distanceNorm;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.valueNames.clear();
    desc.valueNames.push_back("None");
    desc.valueNames.push_back("Sum of Frames");
    desc.valueNames.push_back("Log Sum of Frames");
    list.push_back(desc);
    desc.valueNames.clear();

    desc.identifier = "usespecdiff";
    desc.name = "Use Spectral Difference";
    desc.description = "Whether to use half-wave rectified spectral difference instead of straight spectrum";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = m_defaultFeParams.useSpectralDifference ? 1 : 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "usechroma";
    desc.name = "Use Chroma Frequency Map";
    desc.description = "Whether to use a chroma frequency map instead of the default warped spectrogram";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = m_defaultFeParams.useChromaFrequencyMap ? 1 : 0;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "gradientlimit";
    desc.name = "Gradient Limit";
    desc.description = "Limit of number of frames that will be accepted from one source without a frame from the other source being accepted";
    desc.minValue = 1;
    desc.maxValue = 10;
    desc.defaultValue = m_defaultParams.maxRunCount;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    list.push_back(desc);

    desc.identifier = "zonewidth";
    desc.name = "Search Zone Width";
    desc.description = "Width of the search zone (error margin) either side of the ongoing match position, in seconds";
    desc.minValue = 1;
    desc.maxValue = 60;
    desc.defaultValue = (float)m_defaultParams.blockTime;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.unit = "s";
    list.push_back(desc);

    desc.identifier = "smooth";
    desc.name = "Smooth Path";
    desc.description = "Smooth the path by replacing steps with diagonals";
    desc.minValue = 0;
    desc.maxValue = 1;
    desc.defaultValue = 1;
    desc.isQuantized = true;
    desc.quantizeStep = 1;
    desc.unit = "";
    list.push_back(desc);

    return list;
}

float
MatchVampPlugin::getParameter(std::string name) const
{
    if (name == "serialise") {
        return m_serialise ? 1.0 : 0.0; 
    } else if (name == "framenorm") {
        return (int)m_feParams.frameNorm;
    } else if (name == "distnorm") {
        return (int)m_params.distanceNorm;
    } else if (name == "usespecdiff") {
        return m_feParams.useSpectralDifference ? 1.0 : 0.0;
    } else if (name == "usechroma") {
        return m_feParams.useChromaFrequencyMap ? 1.0 : 0.0;
    } else if (name == "gradientlimit") {
        return m_params.maxRunCount;
    } else if (name == "zonewidth") {
        return (float)m_params.blockTime;
    } else if (name == "smooth") {
        return m_smooth ? 1.0 : 0.0;
    }
    
    return 0.0;
}

void
MatchVampPlugin::setParameter(std::string name, float value)
{
    if (name == "serialise") {
        m_serialise = (value > 0.5);
    } else if (name == "framenorm") {
        m_feParams.frameNorm = (FeatureExtractor::FrameNormalisation)(int(value + 0.1));
    } else if (name == "distnorm") {
        m_params.distanceNorm = (DistanceMetric::DistanceNormalisation)(int(value + 0.1));
    } else if (name == "usespecdiff") {
        m_feParams.useSpectralDifference = (value > 0.5);
    } else if (name == "usechroma") {
        m_feParams.useChromaFrequencyMap = (value > 0.5);
    } else if (name == "gradientlimit") {
        m_params.maxRunCount = int(value + 0.1);
    } else if (name == "zonewidth") {
        m_params.blockTime = value;
    } else if (name == "smooth") {
        m_smooth = (value > 0.5);
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
    m_params.fftSize = m_blockSize;
    m_feParams.fftSize = m_blockSize;
    m_fe1 = new FeatureExtractor(m_feParams);
    m_fe2 = new FeatureExtractor(m_feParams);
    m_pm1 = new Matcher(m_params, 0, m_fe1->getFeatureSize());
    m_pm2 = new Matcher(m_params, m_pm1, m_fe2->getFeatureSize());
    m_pm1->setOtherMatcher(m_pm2);
    m_feeder = new MatchFeatureFeeder(m_pm1, m_pm2);
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

    m_stepSize = stepSize;
    m_stepTime = float(stepSize) / m_inputSampleRate;
    m_blockSize = blockSize;

    createMatchers();
    m_begin = true;
    m_locked = false;

    return true;
}

void
MatchVampPlugin::reset()
{
    delete m_feeder;
    delete m_fe1;
    delete m_fe2;
    delete m_pm1;
    delete m_pm2;

    m_feeder = 0;
    m_fe1 = 0;
    m_fe2 = 0;
    m_pm1 = 0;
    m_pm2 = 0;

    m_frameNo = 0;
    m_lastFrameIn1 = 0;
    m_lastFrameIn2 = 0;

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
    m_pathOutNo = list.size();
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
    m_abOutNo = list.size();
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
    m_baOutNo = list.size();
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
    m_abDivOutNo = list.size();
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
    m_abRatioOutNo = list.size();
    list.push_back(desc);

    int featureSize = FeatureExtractor(m_feParams).getFeatureSize();
    
    desc.identifier = "a_features";
    desc.name = "A Features";
    desc.description = "Spectral features extracted from performance A";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = featureSize;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_aFeaturesOutNo = list.size();
    list.push_back(desc);

    desc.identifier = "b_features";
    desc.name = "B Features";
    desc.description = "Spectral features extracted from performance B";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = featureSize;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_bFeaturesOutNo = list.size();
    list.push_back(desc);

    return list;
}

bool
MatchVampPlugin::aboveThreshold(const float *frame)
{
    float threshold = 1e-5f;
    float rms = 0.f;
    for (int i = 0; i < m_blockSize/2 + 2; ++i) {
        rms += frame[i] * frame[i];
    }
    rms = sqrtf(rms / (m_blockSize/2 + 2));
    return (rms > threshold);
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

    if (aboveThreshold(inputBuffers[0])) m_lastFrameIn1 = m_frameNo;
    if (aboveThreshold(inputBuffers[1])) m_lastFrameIn2 = m_frameNo;

    vector<double> f1 = m_fe1->process(inputBuffers[0]);
    vector<double> f2 = m_fe2->process(inputBuffers[1]);
    
    m_feeder->feed(f1, f2);

    FeatureSet returnFeatures;

    Feature f;
    f.hasTimestamp = false;

    f.values.clear();
    for (int j = 0; j < (int)f1.size(); ++j) {
        f.values.push_back(float(f1[j]));
    }
    returnFeatures[m_aFeaturesOutNo].push_back(f);

    f.values.clear();
    for (int j = 0; j < (int)f2.size(); ++j) {
        f.values.push_back(float(f2[j]));
    }
    returnFeatures[m_bFeaturesOutNo].push_back(f);

//    std::cerr << ".";
//    std::cerr << std::endl;

    ++m_frameNo;
    
    return returnFeatures;
}

MatchVampPlugin::FeatureSet
MatchVampPlugin::getRemainingFeatures()
{
    m_feeder->finish();

    FeatureSet returnFeatures;
    
    Finder *finder = m_feeder->getFinder();
    std::vector<int> pathx;
    std::vector<int> pathy;
    int len = finder->retrievePath(m_smooth, pathx, pathy);

    int prevx = 0;
    int prevy = 0;

    for (int i = 0; i < len; ++i) {

        int x = pathx[i];
        int y = pathy[i];

        Vamp::RealTime xt = Vamp::RealTime::frame2RealTime
            (x * m_stepSize, lrintf(m_inputSampleRate));
        Vamp::RealTime yt = Vamp::RealTime::frame2RealTime
            (y * m_stepSize, lrintf(m_inputSampleRate));

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

    delete m_feeder;
    delete m_pm1;
    delete m_pm2;
    m_feeder = 0;
    m_pm1 = 0;
    m_pm2 = 0;

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
