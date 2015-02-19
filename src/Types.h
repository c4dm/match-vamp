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

#ifndef MATCH_TYPES_H
#define MATCH_TYPES_H

#include <vector>

/// A single value in a feature vector
typedef double featurebin_t;

/// A feature vector
typedef std::vector<featurebin_t> feature_t;

/// A sequence of feature vectors
typedef std::vector<feature_t> featureseq_t;

/// The distance between two feature vectors
typedef float distance_t;

/// A distance vector
typedef std::vector<distance_t> distancevec_t;

/// A distance matrix
typedef std::vector<distancevec_t> distancemat_t;

/// The integrated distance (path cost) from the origin to a given point
typedef double pathcost_t;

/// A vector of path costs
typedef std::vector<pathcost_t> pathcostvec_t;

/// A matrix of path costs
typedef std::vector<pathcostvec_t> pathcostmat_t;

/// A direction advance instruction or state
enum advance_t {
    AdvanceNone,
    AdvanceBoth,
    AdvanceThis,
    AdvanceOther
};

/// A vector of advance directions
typedef std::vector<advance_t> advancevec_t;

/// A matrix of advance directions
typedef std::vector<advancevec_t> advancemat_t;


#endif
