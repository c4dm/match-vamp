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

using std::vector;

double
DistanceMetric::calcDistance(const vector<double> &f1,
			     const vector<double> &f2)
{
    double d = 0;
    double sum = 0;

    int featureSize = f1.size();
    assert(int(f2.size()) == featureSize);
    
    for (int i = 0; i < featureSize; i++) {
        assert(f1[i] >= 0);
        assert(f2[i] >= 0);
        d += fabs(f1[i] - f2[i]);
        sum += f1[i] + f2[i];
    }

    if (sum == 0)
        return 0;
    if (m_norm == NormaliseDistanceToSum)
        return d / sum; // 0 <= d/sum <= 2
    if (m_norm != NormaliseDistanceToLogSum)
        return d;

    // note if this were to be restored, it would have to use
    // totalEnergies vector instead of f1[freqMapSize] which used to
    // store the total energy:
    //	double weight = (5 + Math.log(f1[freqMapSize] + f2[freqMapSize]))/10.0;

    double weight = (8 + log(sum)) / 10.0;

    if (weight < 0) weight = 0;
    else if (weight > 1) weight = 1;

    return d / sum * weight;
}

