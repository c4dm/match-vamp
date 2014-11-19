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

#ifndef _FINDER_H_
#define _FINDER_H_

#include <vector>
#include <iostream>

#include "Matcher.h"

class Finder
{
public:
    Finder(Matcher *pm);

    ~Finder();

    /**
     * Tell the finder that one or both files ends sooner than it
     * thought, i.e. that some of the trailing features are silence or
     * otherwise to be ignored. d1 and d2 are feature frame counts for
     * matchers 1 and 2 respectively. If this is not called, the full
     * duration of each input will be considered.
     */
    void setDurations(int d1, int d2);
    
    Matcher::Advance getExpandDirection(int row, int col);
    Matcher::Advance getExpandDirection(int row, int col, bool check);
    
    /** Calculates a rectangle of the path cost matrix so that the
     *  minimum cost path between the bottom left and top right
     *  corners can be computed.  Caches previous values to avoid
     *  calling find() multiple times, and is several times faster as
     *  a result.
     *
     *  @param r1 the bottom of the rectangle to be calculated
     *  @param c1 the left side of the rectangle to be calculated
     *  @param r2 the top of the rectangle to be calculated
     *  @param c2 the right side of the rectangle to be calculated
     */
    void recalculatePathCostMatrix(int r1, int c1, int r2, int c2);

    /**
     * Track back after all of the matchers have been fed in order to
     * obtain the lowest cost path available. Path x and y coordinate
     * pairs are returned in corresponding elements of pathx and
     * pathy. Return value is the length of the returned path: only
     * this many elements from pathx and pathy are valid (any
     * subsequent ones may be spurious).
     *
     * @param smooth whether to smooth the path before returning it
     */
    int retrievePath(bool smooth, std::vector<int> &pathx, std::vector<int> &pathy);

protected:
    Matcher *m_m;
    int m_duration1;
    int m_duration2;
};

#endif
