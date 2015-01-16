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

#ifndef DISTANCE_METRIC_H
#define DISTANCE_METRIC_H

#include <vector>

class DistanceMetric
{
public:
    enum DistanceNormalisation {
            
        /** Do not normalise distance metrics */
        NoDistanceNormalisation,

        /** Normalise distance metric for pairs of frames by the sum
         *  of the two frames. */
        NormaliseDistanceToSum,

        /** Normalise distance metric for pairs of frames by the log
         *  of the sum of the frames. */
        NormaliseDistanceToLogSum,
    };

    struct Parameters {

        Parameters() :
            norm(NormaliseDistanceToLogSum)
        {}

        /** Normalisation for distance metrics. */
        DistanceNormalisation norm;
    };
    
    DistanceMetric(Parameters params);
    
    /** Calculates the Manhattan distance between two vectors, with an
     *  optional normalisation by the combined values in the
     *  vectors. Since the vectors contain energy, this could be
     *  considered as a squared Euclidean distance metric. Note that
     *  normalisation assumes the values are all non-negative.
     *
     *  @param f1 one of the vectors involved in the distance calculation
     *  @param f2 one of the vectors involved in the distance calculation
     *  @return the distance
     */
    double calcDistance(const std::vector<double> &f1,
			const std::vector<double> &f2);
    
private:
    Parameters m_params;
};

#endif
