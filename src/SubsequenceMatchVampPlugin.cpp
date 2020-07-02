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

#include "SubsequenceMatchVampPlugin.h"

#include <vamp/vamp.h>
#include <vamp-sdk/RealTime.h>

#include <vector>
#include <algorithm>

using std::string;

// We want to ensure our freq map / crossover bin are always valid
// with a fixed FFT length in seconds, so must reject low sample rates
static float sampleRateMin = 5000.f;

static float defaultStepTime = 0.020f;

SubsequenceMatchVampPlugin::SubsequenceMatchVampPlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(int(inputSampleRate * defaultStepTime + 0.001)),
    m_stepTime(defaultStepTime),
    m_blockSize(2048),
    m_coarseDownsample(50),
    m_serialise(false),
    m_smooth(false),
    m_params(defaultStepTime),
    m_defaultParams(defaultStepTime),
    m_feParams(inputSampleRate),
    m_defaultFeParams(44100), // parameter descriptors can't depend on samplerate
    m_secondReferenceFrequency(m_defaultFeParams.referenceFrequency),
    m_fcParams(),
    m_defaultFcParams(),
    m_dParams(),
    m_defaultDParams()
{
    if (inputSampleRate < sampleRateMin) {
        std::cerr << "SubsequenceMatchVampPlugin::SubsequenceMatchVampPlugin: input sample rate "
                  << inputSampleRate << " < min supported rate "
                  << sampleRateMin << ", plugin will refuse to initialise"
                  << std::endl;
    }
}

SubsequenceMatchVampPlugin::~SubsequenceMatchVampPlugin()
{
}

string
SubsequenceMatchVampPlugin::getIdentifier() const
{
    return "match_subsequence";
}

string
SubsequenceMatchVampPlugin::getName() const
{
    return "Match Subsequence Aligner";
}

string
SubsequenceMatchVampPlugin::getDescription() const
{
    return "Calculate alignment between a reference performance and a performance known to represent only part of the same material";
}

string
SubsequenceMatchVampPlugin::getMaker() const
{
    return "Simon Dixon (plugin by Chris Cannam)";
}

int
SubsequenceMatchVampPlugin::getPluginVersion() const
{
    return 3;
}

string
SubsequenceMatchVampPlugin::getCopyright() const
{
    return "GPL";
}

SubsequenceMatchVampPlugin::ParameterList
SubsequenceMatchVampPlugin::getParameterDescriptors() const
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

#ifdef USE_COMPACT_TYPES
    desc.identifier = "scale";
    desc.name = "Distance scale";
    desc.description = "Scale factor to use when mapping distance metric into byte range for storage";
    desc.minValue = 1;
    desc.maxValue = 1000;
    desc.defaultValue = float(m_defaultDParams.scale);
    desc.isQuantized = false;
    list.push_back(desc);
#endif
    
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

    //!!! + m_coarseDownsample, + reconsider params not useful in this plugin
    
    
    return list;
}

float
SubsequenceMatchVampPlugin::getParameter(std::string name) const
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
    } else if (name == "scale") {
        return float(m_dParams.scale);
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
SubsequenceMatchVampPlugin::setParameter(std::string name, float value)
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
    } else if (name == "scale") {
        m_dParams.scale = value;
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
SubsequenceMatchVampPlugin::getPreferredStepSize() const
{
    return int(m_inputSampleRate * defaultStepTime + 0.001);
}

size_t
SubsequenceMatchVampPlugin::getPreferredBlockSize() const
{
    return m_defaultFeParams.fftSize;
}

bool
SubsequenceMatchVampPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (m_inputSampleRate < sampleRateMin) {
        std::cerr << "SubsequenceMatchVampPlugin::SubsequenceMatchVampPlugin: input sample rate "
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

    m_params.hopTime = m_stepTime;
    m_feParams.fftSize = m_blockSize;

    return true;
}

void
SubsequenceMatchVampPlugin::reset()
{
    //!!!
}

SubsequenceMatchVampPlugin::OutputList
SubsequenceMatchVampPlugin::getOutputDescriptors() const
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
    
    return list;
}

SubsequenceMatchVampPlugin::FeatureSet
SubsequenceMatchVampPlugin::process(const float *const *inputBuffers,
                                    Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;
    
    return returnFeatures;
}

SubsequenceMatchVampPlugin::FeatureSet
SubsequenceMatchVampPlugin::getRemainingFeatures()
{
#ifdef _WIN32
    HANDLE mutex;
#else
    pthread_mutex_t mutex;
#endif

    if (m_serialise) {
#ifdef _WIN32
        WaitForSingleObject(mutex, INFINITE);
#else
        pthread_mutex_lock(&mutex);
#endif
    }

    FeatureSet returnFeatures;

    if (m_serialise) {
#ifdef _WIN32
        ReleaseMutex(mutex);
#else
        pthread_mutex_unlock(&mutex);
#endif
    }

    return returnFeatures;
}
