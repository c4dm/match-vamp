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
#include <map>

using namespace std;

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
    m_params(inputSampleRate, defaultStepTime, m_blockSize),
    m_defaultParams(inputSampleRate, defaultStepTime, m_blockSize),
    m_feParams(inputSampleRate, m_blockSize),
    m_defaultFeParams(inputSampleRate, m_blockSize),
    m_fcParams(),
    m_defaultFcParams()
{
    if (inputSampleRate < sampleRateMin) {
        cerr << "MatchVampPlugin::MatchVampPlugin: input sample rate "
                  << inputSampleRate << " < min supported rate "
                  << sampleRateMin << ", plugin will refuse to initialise"
                  << endl;
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
}

MatchVampPlugin::~MatchVampPlugin()
{
//    cerr << "MatchVampPlugin::~MatchVampPlugin(" << this << "): extant = " << --extant << endl;

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
    desc.defaultValue = (int)m_defaultFcParams.norm;
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
    desc.defaultValue = (int)m_defaultFcParams.order;
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

    desc.identifier = "diagonalweight";
    desc.name = "Diagonal Weight";
    desc.description = "Weight applied to cost of diagonal step relative to horizontal or vertical step. The default of 2.0 is good for gross tracking of quite different performances; closer to 1.0 produces a smoother path for performances more similar in tempo";
    desc.minValue = 1.0;
    desc.maxValue = 2.0;
    desc.defaultValue = 2.0;
    desc.isQuantized = false;
    desc.unit = "";
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
MatchVampPlugin::getParameter(string name) const
{
    if (name == "serialise") {
        return m_serialise ? 1.0 : 0.0; 
    } else if (name == "framenorm") {
        return (int)m_fcParams.norm;
    } else if (name == "distnorm") {
        return (int)m_params.distanceNorm;
    } else if (name == "usespecdiff") {
        return (int)m_fcParams.order;
    } else if (name == "usechroma") {
        return m_feParams.useChromaFrequencyMap ? 1.0 : 0.0;
    } else if (name == "gradientlimit") {
        return m_params.maxRunCount;
    } else if (name == "diagonalweight") {
        return m_params.diagonalWeight;
    } else if (name == "zonewidth") {
        return (float)m_params.blockTime;
    } else if (name == "smooth") {
        return m_smooth ? 1.0 : 0.0;
    }
    
    return 0.0;
}

void
MatchVampPlugin::setParameter(string name, float value)
{
    if (name == "serialise") {
        m_serialise = (value > 0.5);
    } else if (name == "framenorm") {
        m_fcParams.norm = (FeatureConditioner::Normalisation)(int(value + 0.1));
    } else if (name == "distnorm") {
        m_params.distanceNorm = (DistanceMetric::DistanceNormalisation)(int(value + 0.1));
    } else if (name == "usespecdiff") {
        m_fcParams.order = (FeatureConditioner::OutputOrder)(int(value + 0.1));
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

    m_pipeline = new MatchPipeline(m_feParams, m_fcParams, m_params);
}

bool
MatchVampPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_inputSampleRate < sampleRateMin) {
        cerr << "MatchVampPlugin::MatchVampPlugin: input sample rate "
                  << m_inputSampleRate << " < min supported rate "
                  << sampleRateMin << endl;
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
    delete m_pipeline;
    m_pipeline = 0;
    m_frameNo = 0;

    m_mag1.clear();
    m_mag2.clear();
    
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

    desc.identifier = "feature_distance";
    desc.name = "Feature Distance";
    desc.description = "Value of distance metric between features at each point-in-A along the chosen alignment path";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::VariableSampleRate;
    desc.sampleRate = outRate;
    m_distOutNo = list.size();
    list.push_back(desc);

    desc.identifier = "feature_mag";
    desc.name = "Feature Magnitudes";
    desc.description = "Sum of magnitudes of feature pairs for each point-in-A along the chosen alignment path";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_magOutNo = list.size();
    list.push_back(desc);

    desc.identifier = "confidence";
    desc.name = "Confidence";
    desc.description = "Confidence metric for the quality of match at each point-in-A along the chosen alignment path";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 1;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_confidenceOutNo = list.size();
    list.push_back(desc);

    desc.identifier = "confidence_peaks";
    desc.name = "Confidence Peaks";
    desc.description = "Peak locations for the confidence metric";
    desc.unit = "";
    desc.hasFixedBinCount = true;
    desc.binCount = 0;
    desc.hasKnownExtents = false;
    desc.isQuantized = false;
    desc.sampleType = OutputDescriptor::FixedSampleRate;
    desc.sampleRate = outRate;
    m_confPeakOutNo = list.size();
    list.push_back(desc);

    return list;
}

static double
magOf(const vector<double> &f)
{
    double mag = 0.0;
    for (int j = 0; j < (int)f.size(); ++j) {
        mag += f[j] * f[j];
    }
    mag = sqrt(mag);
    return mag;
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
    
//    cerr << timestamp.toString();

    m_pipeline->feedFrequencyDomainAudio(inputBuffers[0], inputBuffers[1]);

    vector<double> f1, f2, c1, c2;

    m_pipeline->extractFeatures(f1, f2);
    m_pipeline->extractConditionedFeatures(c1, c2);    

    m_mag1.push_back(magOf(f1));
    m_mag2.push_back(magOf(f2));

    FeatureSet returnFeatures;

    Feature f;
    f.hasTimestamp = false;

    f.values.clear();
    for (int j = 0; j < (int)c1.size(); ++j) {
        f.values.push_back(float(c1[j]));
    }
    returnFeatures[m_aFeaturesOutNo].push_back(f);

    f.values.clear();
    for (int j = 0; j < (int)c2.size(); ++j) {
        f.values.push_back(float(c2[j]));
    }
    returnFeatures[m_bFeaturesOutNo].push_back(f);

//    cerr << ".";
//    cerr << endl;

    ++m_frameNo;
    
    return returnFeatures;
}

MatchVampPlugin::FeatureSet
MatchVampPlugin::getRemainingFeatures()
{
    m_pipeline->finish();

    FeatureSet returnFeatures;
    
    Finder *finder = m_pipeline->getFinder();
    vector<int> pathx;
    vector<int> pathy;
    vector<float> distances;
    finder->retrievePath(pathx, pathy, distances);

    int prevx = 0;
    int prevy = 0;
    int len = pathx.size();

//!!!
//  m_smooth = true;
    
    if (m_smooth) {
    
        vector<float> confidence;

        for (int i = 0; i < len; ++i) {
            int x = pathx[i];
            int y = pathy[i];
            if (x != prevx) {
                double magSum = m_mag1[y] + m_mag2[x];
                double distance = distances[i];
                float c = magSum - distance;
                confidence.push_back(c);

                Feature f;
                f.values.push_back(c);
                returnFeatures[m_confidenceOutNo].push_back(f);

                f.values.clear();
                f.values.push_back(magSum);
                returnFeatures[m_magOutNo].push_back(f);

                f.values.clear();
                f.values.push_back(distance);
                returnFeatures[m_distOutNo].push_back(f);
            }
        }
        
        if (!confidence.empty()) {

            map<int, int> pinpoints;
            
            vector<float> csorted = confidence;
            sort(csorted.begin(), csorted.end());
            float thresh = csorted[int(csorted.size() * 0.7)]; // 70th percentile

            for (int i = 1; i + 1 < int(confidence.size()); ++i) {

                int x = pathx[i];
                int y = pathy[i];

                if (confidence[i] > thresh &&
                    confidence[i] > confidence[i-1] &&
                    confidence[i] >= confidence[i+1]) {

                    pinpoints[x] = y;

                    Vamp::RealTime xt = Vamp::RealTime::frame2RealTime
                        (x * m_stepSize, lrintf(m_inputSampleRate));
                    Feature feature;
                    feature.hasTimestamp = true;
                    feature.timestamp = m_startTime + xt;
                    returnFeatures[m_confPeakOutNo].push_back(feature);
                }
            }

            finder->smoothWithPinPoints(pinpoints);
        }

        pathx.clear();
        pathy.clear();
        distances.clear();
        finder->retrievePath(pathx, pathy, distances);
        len = pathx.size();
    }    

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
        cerr << i << ": [" << pathx[i] << "," << pathy[i] << "]" << endl;
    }

    cerr << endl;
    cerr << "File: A" << endl;
    cerr << "Marks: -1" << endl;
    cerr << "FixedPoints: true 0" << endl;
    cerr << "0" << endl;
    cerr << "0" << endl;
    cerr << "0" << endl;
    cerr << "0" << endl;
    cerr << "File: B" << endl;
    cerr << "Marks: 0" << endl;
    cerr << "FixedPoints: true 0" << endl;
    cerr << "0.02" << endl;
    cerr << "0.02" << endl;

    cerr << len << endl;
    for (int i = 0; i < len; ++i) {
        cerr << pathx[i] << endl;
    }

    cerr << len << endl;
    for (int i = 0; i < len; ++i) {
        cerr << pathy[i] << endl;
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
