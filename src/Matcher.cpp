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

#include "Matcher.h"

#include <iostream>

#include <cstdlib>
#include <cassert>

using namespace std;

//#define DEBUG_MATCHER 1

Matcher::Matcher(Parameters parameters, Matcher *p, int m_featureSize_) :
    m_params(parameters),
    m_featureSize(m_featureSize_),
    m_metric(parameters.distanceNorm)
{
#ifdef DEBUG_MATCHER
    cerr << "Matcher::Matcher(" << m_params.sampleRate << ", " << p << ", " << m_featureSize << ")" << endl;
#endif

    m_otherMatcher = p;	// the first matcher will need this to be set later
    m_firstPM = (!p);
    m_frameCount = 0;
    m_runCount = 0;
    m_blockSize = 0;

    m_blockSize = lrint(m_params.blockTime / m_params.hopTime);
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

    m_frames = vector<vector<double> >
        (m_blockSize, vector<double>(m_featureSize, -1.0));

    m_distXSize = m_blockSize * 2;

    size();

    m_frameCount = 0;
    m_runCount = 0;
    
    m_initialised = true;
}

void
Matcher::size()
{
    int distSize = (m_params.maxRunCount + 1) * m_blockSize;
    m_bestPathCost.resize(m_distXSize, vector<double>(distSize, -1));
    m_distance.resize(m_distXSize, vector<float>(distSize, -1));
    m_advance.resize(m_distXSize, vector<Advance>(distSize, AdvanceNone));
    m_first.resize(m_distXSize, 0);
    m_last.resize(m_distXSize, 0);
}

void
Matcher::consumeFeatureVector(vector<double> feature)
{
    if (!m_initialised) init();
    int frameIndex = m_frameCount % m_blockSize; 
    m_frames[frameIndex] = feature;
    calcAdvance();
}

void
Matcher::calcAdvance()
{
    int frameIndex = m_frameCount % m_blockSize;

    if (m_frameCount >= m_distXSize) {
        m_distXSize *= 2;
        size();
    }

    if (m_firstPM && (m_frameCount >= m_blockSize)) {

        int len = m_last[m_frameCount - m_blockSize] -
                 m_first[m_frameCount - m_blockSize];

        // We need to copy distance[m_frameCount-m_blockSize] to
        // distance[m_frameCount], and then truncate
        // distance[m_frameCount-m_blockSize] to its first len elements.
        // Same for bestPathCost.

        vector<float> dOld = m_distance[m_frameCount - m_blockSize];
        vector<float> dNew(len, -1.f);

        vector<double> bpcOld = m_bestPathCost[m_frameCount - m_blockSize];
        vector<double> bpcNew(len, -1.0);

        vector<Advance> adOld = m_advance[m_frameCount - m_blockSize];
        vector<Advance> adNew(len, AdvanceNone);

        for (int i = 0; i < len; ++i) {
            dNew[i] = dOld[i];
            bpcNew[i] = bpcOld[i];
            adNew[i] = adOld[i];
        }
        
        m_distance[m_frameCount] = dOld;
        m_distance[m_frameCount - m_blockSize] = dNew;

        m_bestPathCost[m_frameCount] = bpcOld;
        m_bestPathCost[m_frameCount - m_blockSize] = bpcNew;

        m_advance[m_frameCount] = adOld;
        m_advance[m_frameCount - m_blockSize] = adNew;
    }

    int stop = m_otherMatcher->m_frameCount;
    int index = stop - m_blockSize;
    if (index < 0)
        index = 0;
    m_first[m_frameCount] = index;
    m_last[m_frameCount] = stop;

    float mn= -1;
    float mx= -1;
    for ( ; index < stop; index++) {

        float dMN = (float) m_metric.calcDistance
            (m_frames[frameIndex],
             m_otherMatcher->m_frames[index % m_blockSize]);
        
        if (mx<0)
            mx = mn = dMN;
        else if (dMN > mx)
            mx = dMN;
        else if (dMN < mn)
            mn = dMN;

        if ((m_frameCount == 0) && (index == 0))    // first element
            updateValue(0, 0, AdvanceNone, 0, dMN);
        else if (m_frameCount == 0)                 // first row
            updateValue(0, index, AdvanceOther,
                        getPathCost(0, index-1), dMN);
        else if (index == 0)                      // first column
            updateValue(m_frameCount, index, AdvanceThis,
                        getPathCost(m_frameCount - 1, 0), dMN);
        else if (index == m_otherMatcher->m_frameCount - m_blockSize) {
            // missing value(s) due to cutoff
            //  - no previous value in current row (resp. column)
            //  - no diagonal value if prev. dir. == curr. dirn
            double min2 = getPathCost(m_frameCount - 1, index);
            //	if ((m_firstPM && (first[m_frameCount - 1] == index)) ||
            //			(!m_firstPM && (m_last[index-1] < m_frameCount)))
            if (m_first[m_frameCount - 1] == index)
                updateValue(m_frameCount, index, AdvanceThis, min2, dMN);
            else {
                double min1 = getPathCost(m_frameCount - 1, index - 1);
                if (min1 + dMN <= min2)
                    updateValue(m_frameCount, index, AdvanceBoth, min1,dMN);
                else
                    updateValue(m_frameCount, index, AdvanceThis, min2,dMN);
            }
        } else {
            double min1 = getPathCost(m_frameCount, index-1);
            double min2 = getPathCost(m_frameCount - 1, index);
            double min3 = getPathCost(m_frameCount - 1, index-1);
            if (min1 <= min2) {
                if (min3 + dMN <= min1)
                    updateValue(m_frameCount, index, AdvanceBoth, min3,dMN);
                else
                    updateValue(m_frameCount, index, AdvanceOther,min1,dMN);
            } else {
                if (min3 + dMN <= min2)
                    updateValue(m_frameCount, index, AdvanceBoth, min3,dMN);
                else
                    updateValue(m_frameCount, index, AdvanceThis, min2,dMN);
            }
        }
        m_otherMatcher->m_last[index]++;
    } // loop for row (resp. column)

    m_frameCount++;
    m_runCount++;

    m_otherMatcher->m_runCount = 0;
}

bool
Matcher::isInRange(int i, int j)
{
    if (m_firstPM) {
        return ((i >= 0) &&
                (i < int(m_first.size())) &&
                (j >= m_first[i]) &&
                (j < int(m_first[i] + m_bestPathCost[i].size())));
    } else {
        return m_otherMatcher->isInRange(j, i);
    }
}

bool
Matcher::isAvailable(int i, int j)
{
    if (m_firstPM) {
        if (isInRange(i, j)) {
            return (m_bestPathCost[i][j - m_first[i]] >= 0);
        } else {
            return false;
        }
    } else {
        return m_otherMatcher->isAvailable(j, i);
    }
}

pair<int, int>
Matcher::getColRange(int i)
{
    if (i < 0 || i >= int(m_first.size())) {
        cerr << "ERROR: Matcher::getColRange(" << i << "): Index out of range"
             << endl;
        throw "Index out of range";
    } else {
        return pair<int, int>(m_first[i], m_last[i]);
    }
}

pair<int, int>
Matcher::getRowRange(int i)
{
    return m_otherMatcher->getColRange(i);
}

float
Matcher::getDistance(int i, int j)
{
    if (m_firstPM) {
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::getDistance(" << i << ", " << j << "): "
                 << "Location is not in range" << endl;
            throw "Distance not available";
        }
        float dist = m_distance[i][j - m_first[i]];
        if (dist < 0) {
            cerr << "ERROR: Matcher::getDistance(" << i << ", " << j << "): "
                 << "Location is in range, but distance ("
                 << dist << ") is invalid or has not been set" << endl;
            throw "Distance not available";
        }
        return dist;
    } else {
        return m_otherMatcher->getDistance(j, i);
    }
}
                
void
Matcher::setDistance(int i, int j, float distance)
{
    if (m_firstPM) {
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::setDistance(" << i << ", " << j << ", "
                 << distance << "): Location is out of range" << endl;
            throw "Indices out of range";
        }
        m_distance[i][j - m_first[i]] = distance;
    } else {
        m_otherMatcher->setDistance(j, i, distance);
    }
}

double
Matcher::getPathCost(int i, int j)
{
    if (m_firstPM) {
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
        return m_bestPathCost[i][j - m_first[i]];
    } else {
        return m_otherMatcher->getPathCost(j, i);
    }
}
                
void
Matcher::setPathCost(int i, int j, Advance dir, double pathCost)
{
    if (m_firstPM) {
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::setPathCost(" << i << ", " << j << ", "
                 << dir << ", " << pathCost
                 << "): Location is out of range" << endl;
            throw "Indices out of range";
        }
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
Matcher::updateValue(int i, int j, Advance dir, double value, float dMN)
{
    if (m_firstPM) {

        m_distance[i][j - m_first[i]] = dMN;
        setPathCost(i, j, dir, value + (dir == AdvanceBoth ? dMN*2: dMN));

    } else {

        int idx = i - m_otherMatcher->m_first[j];
        
        if (idx == (int)m_otherMatcher->m_distance[j].size()) {
            // This should never happen, but if we allow arbitrary
            // pauses in either direction, and arbitrary lengths at
            // end, it is better than a segmentation fault.
            cerr << "Emergency resize: " << idx << " -> " << idx * 2 << endl;
            m_otherMatcher->m_bestPathCost[j].resize(idx * 2, -1);
            m_otherMatcher->m_distance[j].resize(idx * 2, -1);
            m_otherMatcher->m_advance[j].resize(idx * 2, AdvanceNone);
        }

        m_otherMatcher->m_distance[j][idx] = dMN;
        m_otherMatcher->setPathCost(j, i, dir, value + (dir == AdvanceBoth ? dMN*2: dMN));
    }
}

Matcher::Advance
Matcher::getAdvance(int i, int j)
{
    if (m_firstPM) {
        if (!isInRange(i, j)) {
            cerr << "ERROR: Matcher::getAdvance(" << i << ", " << j << "): "
                 << "Location is not in range" << endl;
            throw "Advance not available";
        }
        return m_advance[i][j - m_first[i]];
    } else {
        return m_otherMatcher->getAdvance(j, i);
    }
}
