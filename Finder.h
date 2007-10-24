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

/** Maps cost matrix coordinates into an efficient
 *  (linear instead of quadratic space) representation.
 *  Stores result of most recent mapping for fast
 *  sequential access.
 */
class Finder {

protected:
    Matcher *pm1, *pm2;
    int index1, index2, bestRow, bestCol;
    int *rowRange;
    int *colRange;

public:
    Finder(Matcher *p1, Matcher *p2);

    ~Finder();

    /** Sets up the instance variables to point to the given
     *  coordinate in the distance matrix.
     *
     *  @param i1 frameNumber in the first Matcher
     *  @param i2 frameNumber in the second Matcher
     *  @return true iff the point (i2,i1) is represented in the distance matrix
     */
    bool find(int i1, int i2);

    /** Returns the range [lo,hi) of legal column indices for the
     *  given row. */
    void getColRange(int row, int *range);

    /** Returns the range [lo,hi) of legal row indices for the given
     *  column. */
    void getRowRange(int col, int *range);

    int getExpandDirection(int row, int col);
    int getExpandDirection(int row, int col, bool check);
	
    unsigned char getDistance(int row, int col);
    void setDistance(int row, int col, unsigned char b);

    int getPathCost(int row, int col);
    int getRawPathCost(int row, int col);
    void setPathCost(int row, int col, int i);

    unsigned char getDistance();
    void setDistance(int b);

    int getPathCost();
    void setPathCost(int i);

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

}; // class Finder

#endif
