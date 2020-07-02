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

#include "DistanceMetric.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

//#define DEBUG_DISTANCE_METRIC 1

template <> uint8_t
DistanceMetric::scaleIntoRange(double distance)
{
    double scaled = m_params.scale * distance;
    if (scaled < 0) {
        scaled = 0;
    }
    if (scaled > DISTANCE_MAX) {
        scaled = DISTANCE_MAX;
        ++m_overcount;
    }
    return uint8_t(scaled);
}

template <> float
DistanceMetric::scaleIntoRange(double distance)
{
    return float(distance);
}

template <> double
DistanceMetric::scaleIntoRange(double distance)
{
    return distance;
}

DistanceMetric::DistanceMetric(Parameters params) :
    m_params(params),
    m_max(0),
    m_overcount(0)
{
#ifdef DEBUG_DISTANCE_METRIC
    cerr << "*** DistanceMetric: metric = " << m_params.metric
         << ", norm = " << m_params.norm
         << ", noise = " << m_params.noise
         << ", scale = " << m_params.scale
         << endl;
#endif
}

DistanceMetric::~DistanceMetric()
{
#ifdef DEBUG_DISTANCE_METRIC
    cerr << "*** DistanceMetric::~DistanceMetric: metric = " << m_params.metric
         << ", norm = " << m_params.norm
         << ", noise = " << m_params.noise;
#ifdef USE_COMPACT_TYPES
    cerr << ", scale = " << m_params.scale;
    cerr << "\n*** DistanceMetric::~DistanceMetric: max scaled value = "
         << distance_print_t(m_max)
         << ", " << m_overcount << " clipped" << endl;
#else
    cerr << ", no scaling";
    cerr << "\n*** DistanceMetric::~DistanceMetric: max value = "
         << distance_print_t(m_max)
         << endl;
#endif
#endif
}

distance_t
DistanceMetric::scaleValueIntoDistanceRange(double value)
{
    return scaleIntoRange<distance_t>(value);
}

distance_t
DistanceMetric::scaleAndTally(double value)
{
    distance_t dist = scaleIntoRange<distance_t>(value);
    if (dist > m_max) m_max = dist;
    return dist;
}

distance_t
DistanceMetric::calcDistance(const feature_t &f1,
			     const feature_t &f2)
{
    double d = 0;
    double sum = 0;
    double eps = 1e-16;

    assert(f2.size() == f1.size());
    int featureSize = static_cast<int>(f1.size());

    double minNoise = 0.0;
#ifdef USE_COMPACT_TYPES
    minNoise = 1.0 / m_params.scale;
#endif
    
    if (m_params.metric == Cosine) {

        double num = 0, denom1 = 0, denom2 = 0;
        
        for (int i = 0; i < featureSize; ++i) {
            num += f1[i] * f2[i];
            denom1 += f1[i] * f1[i];
            denom2 += f2[i] * f2[i];
        }

        d = 1.0 - (num / (eps + sqrt(denom1 * denom2)));

        if (m_params.noise == AddNoise) {
            double noise = 1e-2;
            if (noise < minNoise) noise = minNoise;
            d += noise;
        }
        if (d > 1.0) d = 1.0;
        
        return scaleAndTally(d); // normalisation param ignored
    }

    if (m_params.metric == Manhattan) {
        for (int i = 0; i < featureSize; i++) {
            d += fabs(f1[i] - f2[i]);
            sum += fabs(f1[i]) + fabs(f2[i]);
        }
    } else {
        // Euclidean
        for (int i = 0; i < featureSize; i++) {
            d += (f1[i] - f2[i]) * (f1[i] - f2[i]);
            sum += fabs(f1[i]) + fabs(f2[i]);
        }
        d = sqrt(d);
    }

    if (m_params.noise == AddNoise) {
        double noise = 1e-3 * featureSize;
        if (noise < minNoise) noise = minNoise;
        d += noise;
        sum += noise;
    }
    
    if (sum == 0) {
        return scaleAndTally(0);
    }

    double distance = 0;

    if (m_params.norm == NormaliseDistanceToSum) {

        distance = d / sum; // 0 <= d/sum <= 2

    } else if (m_params.norm == NormaliseDistanceToLogSum) {

        // note if this were to be restored, it would have to use
        // totalEnergies vector instead of f1[freqMapSize] which used to
        // store the total energy:
        //	double weight = (5 + Math.log(f1[freqMapSize] + f2[freqMapSize]))/10.0;

        double weight = (8 + log(sum)) / 10.0;
    
        if (weight < 0) weight = 0;
        else if (weight > 1) weight = 1;

        distance = d / sum * weight;

    } else {

        distance = d;
    }
    
    return scaleAndTally(distance);
}
