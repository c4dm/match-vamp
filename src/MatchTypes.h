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
#include <cstdint>

#ifdef USE_COMPACT_TYPES

/// A single value in a feature vector
typedef float featurebin_t;

/// The distance between two feature vectors
typedef uint8_t distance_t;

const int MaxDistance = 0xfe;
const int InvalidDistance = 0xff;

/// The integrated distance (path cost) from the origin to a given point
typedef uint32_t pathcost_t;

const int MaxPathCost = 0xfffffffe;
const int InvalidPathCost = 0xffffffff;

/// A direction advance instruction or state
enum advance_t : uint8_t {
    AdvanceNone,
    AdvanceBoth,
    AdvanceThis,
    AdvanceOther
};

#else

#ifndef USE_PRECISE_TYPES
#error You must define either USE_COMPACT_TYPES or USE_PRECISE_TYPES
#endif

/// A single value in a feature vector
typedef double featurebin_t;

/// The distance between two feature vectors
typedef float distance_t;

const float MaxDistance = FLOAT_MAX;
const float InvalidDistance = -1.f;

/// The integrated distance (path cost) from the origin to a given point
typedef double pathcost_t;

const double MaxPathCost = DOUBLE_MAX;
const double InvalidPathCost = -1.;

/// A direction advance instruction or state
enum advance_t {
    AdvanceNone,
    AdvanceBoth,
    AdvanceThis,
    AdvanceOther
};

#endif


/// A feature vector
typedef std::vector<featurebin_t> feature_t;

/// A sequence of feature vectors
typedef std::vector<feature_t> featureseq_t;

/// A distance vector
typedef std::vector<distance_t> distancevec_t;

/// A distance matrix
typedef std::vector<distancevec_t> distancemat_t;

/// A vector of path costs
typedef std::vector<pathcost_t> pathcostvec_t;

/// A matrix of path costs
typedef std::vector<pathcostvec_t> pathcostmat_t;

/// A vector of advance directions
typedef std::vector<advance_t> advancevec_t;

/// A matrix of advance directions
typedef std::vector<advancevec_t> advancemat_t;


#endif
