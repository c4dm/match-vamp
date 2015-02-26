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

#include "DistanceMetric.h"

#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

//#define DEBUG_DISTANCE_METRIC 1

template <> uint8_t
DistanceMetric::scaleIntoRange(double distance)
{
    return uint8_t(m_params.scale * distance);
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
    m_params(params)
{
#ifdef DEBUG_DISTANCE_METRIC
    cerr << "*** DistanceMetric: norm = " << m_params.norm
         << endl;
#endif
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

    if (m_params.metric == Cosine) {

        double num = 0, denom1 = 0, denom2 = 0;
        
        for (int i = 0; i < featureSize; ++i) {
            num += f1[i] * f2[i];
            denom1 += f1[i] * f1[i];
            denom2 += f2[i] * f2[i];
        }

        d = 1.0 - (num / (eps + sqrt(denom1 * denom2)));

        if (m_params.noise == AddNoise) {
            d += 1e-2;
        }
        if (d > 1.0) d = 1.0;
        
        return scaleIntoRange<distance_t>(d); // normalisation param ignored
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

    double noise = 1e-3 * featureSize;
    if (m_params.noise == AddNoise) {
        d += noise;
        sum += noise;
    }
    
    if (sum == 0) {
        return scaleIntoRange<distance_t>(0);
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
    
    return scaleIntoRange<distance_t>(distance);
}
