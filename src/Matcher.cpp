/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright (c) 2007-2015 Simon Dixon, Chris Cannam, and Queen Mary
    University of London, Copyright (c) 2014-2015 Tido GmbH.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "Matcher.h"

#include <iostream>

#include <cstdlib>
#include <cassert>

using namespace std;

//#define DEBUG_MATCHER 1
//#define PERFORM_ERROR_CHECKS 1

Matcher::Matcher(Parameters parameters, DistanceMetric::Parameters dparams,
                 Matcher *p) :
    m_params(parameters),
    m_metric(dparams)
{
#ifdef DEBUG_MATCHER
    cerr << "*** Matcher: hopTime = " << parameters.hopTime
         << ", blockTime = " << parameters.blockTime
         << ", maxRunCount = " << parameters.maxRunCount
         << ", diagonalWeight = " << parameters.diagonalWeight << endl;
#endif
    
    m_otherMatcher = p;	// the first matcher will need this to be set later
    m_firstPM = (!p);
    m_frameCount = 0;
    m_runCount = 0;
    m_blockSize = 0;
    m_distXSize = 0;

    m_blockSize = int(m_params.blockTime / m_params.hopTime + 0.5);
#ifdef DEBUG_MATCHER
    cerr << "Matcher: m_blockSize = " << m_blockSize << endl;
#endif

    m_initialised = false;
} 

Matcher::~Matcher()
{
#ifdef DEBUG_MATCHER
    cerr << "Matcher(" << this << ")::~Matcher()" << endl;
#endif
}

void
Matcher::init()
{
    if (m_initialised) return;

    m_features = featureseq_t(m_blockSize);

    m_distXSize = m_blockSize * 2;

    size();

    m_frameCount = 0;
    m_runCount = 0;
    
    m_initialised = true;
}

bool
Matcher::isAvailable(int i, int j)
{
    if (m_firstPM) {
        if (isInRange(i, j)) {
            return (m_distance[i][j - m_first[i]] != INVALID_DISTANCE);
        } else {
            return false;
        }
    } else {
        return m_otherMatcher->isAvailable(j, i);
    }
}

bool
Matcher::isRowAvailable(int i)
{
    if (m_firstPM) {

        if (i < 0 || i >= int(m_first.size())) return false;
        for (auto c: m_distance[i]) {
            if (c != INVALID_DISTANCE) return true;
        }
        return false;

    } else {
        return m_otherMatcher->isColAvailable(i);
    }
}

bool
Matcher::isColAvailable(int j)
{
    if (m_firstPM) {
        for (int i = 0; i < int(m_first.size()); ++i) {
            if (j >= m_first[i] && j < m_last[i]) {
                if (m_distance[i][j - m_first[i]] != INVALID_DISTANCE) {
                    return true;
                }
            }
        }
        return false;
    } else {
        return m_otherMatcher->isRowAvailable(j);
    }
}

bool
Matcher::isInRange(int i, int j)
{
    if (m_firstPM) {
        return ((i >= 0) &&
                (i < int(m_first.size())) &&
                (j >= m_first[i]) &&
                (j < int(m_first[i] + m_distance[i].size())));
    } else {
        return m_otherMatcher->isInRange(j, i);
    }
}

pair<int, int>
Matcher::getColRangeForRow(int i)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (i < 0 || i >= int(m_first.size())) {
            cerr << "ERROR: Matcher::getColRangeForRow(" << i << "): Index out of range"
                 << endl;
            throw "Index out of range";
        }
#endif
        return pair<int, int>(m_first[i], m_last[i]);
    } else {
        return m_otherMatcher->getRowRangeForCol(i);
    }
}

pair<int, int>
Matcher::getRowRangeForCol(int i)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (i < 0 || i >= int(m_otherMatcher->m_first.size())) {
            cerr << "ERROR: Matcher::getRowRangeForCol(" << i << "): Index out of range"
                 << endl;
            throw "Index out of range";
        }
#endif
        return pair<int, int>(m_otherMatcher->m_first[i],
                              m_otherMatcher->m_last[i]);
    } else {
        return m_otherMatcher->getColRangeForRow(i);
    }
}

distance_t
Matcher::getDistance(int i, int j)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::getDistance(" << i << ", " << j << "): "
                 << "Location is not in range" << endl;
            throw "Distance not available";
        }
#endif
        distance_t dist = m_distance[i][j - m_first[i]];
#ifdef PERFORM_ERROR_CHECKS
        if (dist == INVALID_DISTANCE) {
            cerr << "ERROR: Matcher::getDistance(" << i << ", " << j << "): "
                 << "Location is in range, but distance ("
                 << distance_print_t(dist)
                 << ") is invalid or has not been set" << endl;
            throw "Distance not available";
        }
#endif
        return dist;
    } else {
        return m_otherMatcher->getDistance(j, i);
    }
}
                
void
Matcher::setDistance(int i, int j, distance_t distance)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::setDistance(" << i << ", " << j << ", "
                 << distance_print_t(distance)
                 << "): Location is out of range" << endl;
            throw "Indices out of range";
        }
#endif
        m_distance[i][j - m_first[i]] = distance;
    } else {
        m_otherMatcher->setDistance(j, i, distance);
    }
}

normpathcost_t
Matcher::getNormalisedPathCost(int i, int j)
{
    // normalised for path length. 1+ prevents division by zero here
    return normpathcost_t(getPathCost(i, j)) / normpathcost_t(1 + i + j);
}

pathcost_t
Matcher::getPathCost(int i, int j)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (!isAvailable(i, j)) {
            if (!isInRange(i, j)) {
                cerr << "ERROR: Matcher::getPathCost(" << i << ", " << j << "): "
                     << "Location is not in range" << endl;
            } else {
                cerr << "ERROR: Matcher::getPathCost(" << i << ", " << j << "): "
                     << "Location is in range, but pathCost ("
                     << m_bestPathCost[i][j - m_first[i]]
                     << ") is invalid or has not been set" << endl;
            }
            throw "Path cost not available";
        }
#endif
        return m_bestPathCost[i][j - m_first[i]];
    } else {
        return m_otherMatcher->getPathCost(j, i);
    }
}
                
void
Matcher::setPathCost(int i, int j, advance_t dir, pathcost_t pathCost)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::setPathCost(" << i << ", " << j << ", "
                 << dir << ", " << pathCost
                 << "): Location is out of range" << endl;
            throw "Indices out of range";
        }
#endif
        m_advance[i][j - m_first[i]] = dir;
        m_bestPathCost[i][j - m_first[i]] = pathCost;
    } else {
        if (dir == AdvanceThis) {
            dir = AdvanceOther;
        } else if (dir == AdvanceOther) {
            dir = AdvanceThis;
        }
        m_otherMatcher->setPathCost(j, i, dir, pathCost);
    }
}

void
Matcher::size()
{
    m_first.resize(m_distXSize, 0);
    m_last.resize(m_distXSize, 0);

    if (m_firstPM) {
        int distSize = (m_params.maxRunCount + 1) * m_blockSize;
        m_bestPathCost.resize(m_distXSize, pathcostvec_t(distSize, INVALID_PATHCOST));
        m_distance.resize(m_distXSize, distancevec_t(distSize, INVALID_DISTANCE));
        m_advance.resize(m_distXSize, advancevec_t(distSize, AdvanceNone));
    }
}

void
Matcher::consumeFeatureVector(const feature_t &feature)
{
    if (!m_initialised) init();
    int frameIndex = m_frameCount % m_blockSize; 
    m_features[frameIndex] = feature;
    calcAdvance();
}

pathcost_t
Matcher::addToCost(pathcost_t cost, pathcost_t increment)
{
    if (PATHCOST_MAX - increment < cost) {
        return PATHCOST_MAX;
    } else {
        return cost + pathcost_t(increment);
    }
}

void
Matcher::calcAdvance()
{
    int frameIndex = m_frameCount % m_blockSize;

    if (m_frameCount >= m_distXSize) {
        m_distXSize = int(m_distXSize * 1.2);
        size();
    }
    
    if (m_firstPM && (m_frameCount >= m_blockSize)) {
        // Memory reduction for old rows
        int oldidx = m_frameCount - m_blockSize;
        int len = m_last[oldidx] - m_first[oldidx];
        m_distance[oldidx].resize(len);
        m_distance[oldidx].shrink_to_fit();
        m_bestPathCost[oldidx].resize(len);
        m_bestPathCost[oldidx].shrink_to_fit();
        m_advance[oldidx].resize(len);
        m_advance[oldidx].shrink_to_fit();
    }

    int stop = m_otherMatcher->m_frameCount;
    int index = stop - m_blockSize;
    if (index < 0) index = 0;

    m_first[m_frameCount] = index;
    m_last[m_frameCount] = stop;

    for ( ; index < stop; index++) {

        distance_t distance = m_metric.calcDistance
            (m_features[frameIndex],
             m_otherMatcher->m_features[index % m_blockSize]);

        pathcost_t straightIncrement(distance);
        pathcost_t diagIncrement = pathcost_t(distance * m_params.diagonalWeight);

        if ((m_frameCount == 0) && (index == 0)) { // first element

            updateValue(0, 0, AdvanceNone,
                        0,
                        distance);

        } else if (m_frameCount == 0) { // first row

            pathcost_t cost {};

            if (isInRange(0, index - 1)) {
                cost = getPathCost(0, index - 1);
            }
            
            updateValue(0, index, AdvanceOther, cost, distance);
            
        } else if (index == 0) { // first column

            pathcost_t cost {};

            if (isInRange(m_frameCount - 1, 0)) {
                cost = getPathCost(m_frameCount - 1, 0);
            }
            
            updateValue(m_frameCount, index, AdvanceThis, cost, distance);
            
        } else if (index == m_otherMatcher->m_frameCount - m_blockSize) {
            
            // missing value(s) due to cutoff
            //  - no previous value in current row (resp. column)
            //  - no diagonal value if prev. dir. == curr. dirn
            
            pathcost_t min2 = getPathCost(m_frameCount - 1, index);

//            cerr << "NOTE: missing value at i = " << m_frameCount << ", j = "
//                 << index << " (first = " << m_firstPM << ")" << endl;
                
            //	if ((m_firstPM && (first[m_frameCount - 1] == index)) ||
            //			(!m_firstPM && (m_last[index-1] < m_frameCount)))
            if (m_first[m_frameCount - 1] == index) {
                
                updateValue(m_frameCount, index, AdvanceThis,
                            min2, distance);
                
            } else {

                pathcost_t min1 = getPathCost(m_frameCount - 1, index - 1);
                if (addToCost(min1, diagIncrement) <=
                    addToCost(min2, straightIncrement)) {
                    updateValue(m_frameCount, index, AdvanceBoth,
                                min1, distance);
                } else {
                    updateValue(m_frameCount, index, AdvanceThis,
                                min2, distance);
                }
            }

        } else {

            pathcost_t min1 = getPathCost(m_frameCount, index - 1);
            pathcost_t min2 = getPathCost(m_frameCount - 1, index);
            pathcost_t min3 = getPathCost(m_frameCount - 1, index - 1);

            pathcost_t cost1 = addToCost(min1, straightIncrement);
            pathcost_t cost2 = addToCost(min2, straightIncrement);
            pathcost_t cost3 = addToCost(min3, diagIncrement);

            // Choosing is easy if there is a strict cheapest of the
            // three. If two or more share the lowest cost, we choose
            // in order of preference: cost3 (AdvanceBoth), cost2
            // (AdvanceThis), cost1 (AdvanceOther) if we are the first
            // matcher; and cost3 (AdvanceBoth), cost1 (AdvanceOther),
            // cost2 (AdvanceThis) if we are the second matcher.  That
            // is, we always prioritise the diagonal followed by the
            // first matcher.

            if (( m_firstPM && (cost1 <  cost2)) ||
                (!m_firstPM && (cost1 <= cost2))) {
                if (cost3 <= cost1) {
                    updateValue(m_frameCount, index, AdvanceBoth,
                                min3, distance);
                } else {
                    updateValue(m_frameCount, index, AdvanceOther,
                                min1, distance);
                }
            } else {
                if (cost3 <= cost2) {
                    updateValue(m_frameCount, index, AdvanceBoth,
                                min3, distance);
                } else {
                    updateValue(m_frameCount, index, AdvanceThis,
                                min2, distance);
                }
            }
        }
        
        m_otherMatcher->m_last[index]++;
    } // loop for row (resp. column)

    m_frameCount++;
    m_runCount++;

    m_otherMatcher->m_runCount = 0;
}

void
Matcher::updateValue(int i, int j, advance_t dir, pathcost_t value, distance_t distance)
{
    pathcost_t increment = distance;
    if (dir == AdvanceBoth) {
        increment = pathcost_t(increment * m_params.diagonalWeight);
    }

    pathcost_t newValue = addToCost(value, increment);
    if (newValue == PATHCOST_MAX) {
        cerr << "ERROR: Path cost overflow at i=" << i << ", j=" << j << ": "
             << value << " + " << increment << " >= " << PATHCOST_MAX << endl;
        newValue = PATHCOST_MAX;
    }
    
    if (m_firstPM) {

        setDistance(i, j, distance);
        setPathCost(i, j, dir, newValue);

    } else {

        if (dir == AdvanceThis) dir = AdvanceOther;
        else if (dir == AdvanceOther) dir = AdvanceThis;

        int idx = i - m_otherMatcher->m_first[j];
        
        if (idx < 0 || size_t(idx) == m_otherMatcher->m_distance[j].size()) {
            // This should never happen, but if we allow arbitrary
            // pauses in either direction, and arbitrary lengths at
            // end, it is better than a segmentation fault.
//            cerr << "Emergency resize: " << idx << " -> " << idx * 2 << endl;
            m_otherMatcher->m_bestPathCost[j].resize(idx * 2, INVALID_PATHCOST);
            m_otherMatcher->m_distance[j].resize(idx * 2, INVALID_DISTANCE);
            m_otherMatcher->m_advance[j].resize(idx * 2, AdvanceNone);
        }

        m_otherMatcher->setDistance(j, i, distance);
        m_otherMatcher->setPathCost(j, i, dir, newValue);
    }
}

advance_t
Matcher::getAdvance(int i, int j)
{
    if (m_firstPM) {
#ifdef PERFORM_ERROR_CHECKS
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::getAdvance(" << i << ", " << j << "): "
                 << "Location is not in range" << endl;
            throw "Advance not available";
        }
#endif
        return m_advance[i][j - m_first[i]];
    } else {
        return m_otherMatcher->getAdvance(j, i);
    }
}

static double k(size_t sz)
{
    return double(sz) / 1024.0;
}

Matcher::MemoryStats
Matcher::getMemoryStats() const
{
    MemoryStats stats;
    stats.features_k = 0.0;
    if (!m_features.empty()) {
        stats.features_k =
            k(m_features.size() * m_features[0].size() * sizeof(featurebin_t));
    }
    
    size_t cells = 0;
    for (const auto &d: m_distance) {
        cells += d.size();
    }
    
    stats.pathcosts_k = k(cells * sizeof(pathcost_t));
    stats.distances_k = k(cells * sizeof(distance_t));
    stats.advances_k = k(cells * sizeof(advance_t));

    if (m_firstPM && m_otherMatcher) {
        stats = stats + m_otherMatcher->getMemoryStats();
    }
    
    return stats;
}

void
Matcher::printStats()
{
    if (m_firstPM) cerr << endl;
    
    cerr << "Matcher[" << this << "] (" << (m_firstPM ? "first" : "second") << "):" << endl;
    cerr << "- block size " << m_blockSize << ", frame count " << m_frameCount << ", dist x size " << m_distXSize << ", initialised " << m_initialised << endl;

    if (m_features.empty()) {
        cerr << "- have no features yet" << endl;
    } else {
        cerr << "- have " << m_features.size() << " features of " << m_features[0].size() << " bins each" << endl;
    }

    size_t cells = 0;
    for (const auto &d: m_distance) {
        cells += d.size();
    }
    if (m_distance.empty()) {
        cerr << "- have no cells in matrix" << endl;
    } else {
        cerr << "- have " << m_distance.size() << " cols in matrix with avg "
             << double(cells) / double(m_distance.size()) << " rows, total "
             << cells << " cells" << endl;
    }

    if (m_firstPM && m_otherMatcher) {
        m_otherMatcher->printStats();
        MemoryStats stats = getMemoryStats();
        cerr << "Memory: "
             << "features " << stats.features_k << "K, "
             << "path costs " << stats.pathcosts_k << "K, "
             << "distances " << stats.distances_k << "K,\n        "
             << "advances " << stats.advances_k << "K, "
             << "total " << stats.total_k() << "K"
             << endl;
        cerr << endl;
    }
}

