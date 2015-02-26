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

#include "MatchTypes.h"

class DistanceMetric
{
public:
    enum Metric {

        /** Calculate the Manhattan distance between feature
         *  vectors. If the vectors contain energy, as the default
         *  MATCH feature does, this could be considered as a squared
         *  Euclidean distance metric. */
        Manhattan,

        /** Calculate the Euclidean distance between feature vectors. */
        Euclidean,

        /** Calculate the cosine distance between feature vectors. The
         *  normalisation setting will be ignored as the result is
         *  already magnitude-independent. */
        Cosine,
    };

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
    
    enum NoiseAddition {

        /** Don't add noise. */
        NoNoise,

        /** Add a constant noise term. This can help avoid
         *  mis-tracking when one file contains a lot of silence. */
        AddNoise,
    };
    
    struct Parameters {

        Parameters() :
            metric(Manhattan),
            norm(NormaliseDistanceToLogSum),
            noise(AddNoise),
            scale(90.)
        {}

        Metric metric;
        DistanceNormalisation norm;
        NoiseAddition noise;
        double scale;
    };
    
    DistanceMetric(Parameters params);
    
    /** Calculates the distance in some metric between two vectors,
     *  with an optional normalisation by the combined values in the
     *  vectors. Note that normalisation assumes the values are all
     *  non-negative.
     *
     *  @param f1 one of the vectors involved in the distance calculation
     *  @param f2 one of the vectors involved in the distance calculation
     *  @return the distance
     */
    distance_t calcDistance(const feature_t &f1,
                            const feature_t &f2);

    /**
     * Mostly for internal use and testing
     */
    distance_t scaleValueIntoDistanceRange(double value);
    
private:
    Parameters m_params;

    template <typename T> T scaleIntoRange(double);
};

#endif
