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

#ifndef _MATCH_VAMP_PLUGIN_H_
#define _MATCH_VAMP_PLUGIN_H_

#include <vamp-sdk/Plugin.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <pthread.h>
#endif

class Matcher;
class MatchFeeder;

class MatchVampPlugin : public Vamp::Plugin
{
public:
    MatchVampPlugin(float inputSampleRate);
    virtual ~MatchVampPlugin();

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
    void createMatchers() const;
    mutable Matcher *pm1;
    mutable Matcher *pm2;
    mutable MatchFeeder *feeder;
    Vamp::RealTime m_startTime;
    size_t m_stepSize;
    bool m_serialise;
    bool m_begin;
    bool m_locked;

#ifdef _WIN32
    static HANDLE m_serialisingMutex;
#else 
    static pthread_mutex_t m_serialisingMutex;
#endif

    static bool m_serialisingMutexInitialised;
};


#endif
