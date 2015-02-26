/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "FeatureConditioner.h"

#include <iostream>
#include <cmath>

using namespace std;

//#define DEBUG_FEATURE_CONDITIONER 1

FeatureConditioner::FeatureConditioner(Parameters parameters) :
    m_params(parameters),
    m_ltAverage(0.0)
{
#ifdef DEBUG_FEATURE_CONDITIONER
    cerr << "*** FeatureConditioner: norm = " << parameters.norm
         << ", order = " << parameters.order
         << ", silenceThreshold = " << parameters.silenceThreshold
         << ", decay = " << parameters.decay << endl;
#endif
}

feature_t
FeatureConditioner::process(const feature_t &feature)
{
    if (m_prev.empty()) {
	m_prev.resize(feature.size(), 0.0);
    }
    if (m_prev.size() != feature.size()) {
	cerr << "ERROR: FeatureConditioner::process: feature size "
	     << feature.size() << " differs from previous feature size "
	     << m_prev.size() << endl;
	return feature;
    }

    int size = static_cast<int>(feature.size());
    
    feature_t out(size, 0.0);

    double totalEnergy = 0;

    switch (m_params.order) {

    case OutputRectifiedDerivative:
        for (int i = 0; i < size; i++) {
            totalEnergy += feature[i];
            if (feature[i] > m_prev[i]) {
                out[i] = feature[i] - m_prev[i];
            } else {
                out[i] = 0;
            }
        }
        break;

    case OutputDerivative:
        for (int i = 0; i < size; i++) {
            totalEnergy += feature[i];
            out[i] = fabs(feature[i] - m_prev[i]);
        }
        break;
        
    case OutputFeatures:
        for (int i = 0; i < size; i++) {
            totalEnergy += feature[i];
            out[i] = feature[i];
        }
        break;
    }

    if (m_ltAverage == 0.0) {
	m_ltAverage = totalEnergy;
    } else {
	double decay = m_params.decay;
        m_ltAverage = m_ltAverage * decay + totalEnergy * (1.0 - decay);
    }

    if (totalEnergy <= m_params.silenceThreshold) {
        for (int i = 0; i < size; i++) {
            out[i] = 0;
	}
    } else if (m_params.norm == NormaliseToSum1) {
        for (int i = 0; i < size; i++) { 
            out[i] = featurebin_t(out[i] / totalEnergy);
	}
    } else if (m_params.norm == NormaliseToLTAverage) {
        for (int i = 0; i < size; i++) {
            out[i] = featurebin_t(out[i] / m_ltAverage);
	}
    }

    m_prev = feature;
    return out;
}
    
