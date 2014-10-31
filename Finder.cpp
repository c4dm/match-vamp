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

#include "Finder.h"

#include "Path.h"

#include <algorithm>


Finder::Finder(Matcher *p1, Matcher *p2)
{
    if (!p1->firstPM)
        std::cerr << "Warning: wrong args in Finder()" << std::endl;
    pm1 = p1;
    pm2 = p2;
    index1 = 0;
    index2 = 0;
    rowRange = new int[2];
    colRange = new int[2];
} // constructor

Finder::~Finder()
{
    delete[] rowRange;
    delete[] colRange;
}

bool
Finder::find(int i1, int i2)
{
    if (i1 >= 0) {
        index1 = i1;
        index2 = i2 - pm1->first[i1];
    }
    return (i1 >= 0) && (i2 >= pm1->first[i1]) && (i2 < pm1->last[i1]);
} // find()

void
Finder::getColRange(int row, int *range)
{
    range[0] = pm1->first[row];
    range[1] = pm1->last[row];
} // getColRange()

void
Finder::getRowRange(int col, int *range)
{
    range[0] = pm2->first[col];
    range[1] = pm2->last[col];
} // getRowRange()

int
Finder::getExpandDirection(int row, int col)
{
    return getExpandDirection(row, col, false);
} // getExpandDirection()

int
Finder::getExpandDirection(int row, int col, bool check)
{
    int min = getPathCost(row, col);
    bestRow = row;
    bestCol = col;
    getRowRange(col, rowRange);
    if (rowRange[1] > row+1)
        rowRange[1] = row+1;	// don't cheat by looking at future :)
    for (int index = rowRange[0]; index < rowRange[1]; index++) {
        int tmp = getPathCost(index, col);
        if (tmp < min) {
            min = tmp;
            bestRow = index;
        }
    }
    getColRange(row, colRange);
    if (colRange[1] > col+1)
        colRange[1] = col+1;	// don't cheat by looking at future :)
    for (int index = colRange[0]; index < colRange[1]; index++) {
        int tmp = getPathCost(row, index);
        if (tmp < min) {
            min = tmp;
            bestCol = index;
            bestRow = row;
        }
    }
    //	System.err.print("  BEST: " + bestRow + " " + bestCol + " " + check);
    //	System.err.println(" " + pm1->frameCount + " " + pm2->frameCount);
    if (check) {
        //		System.err.println(find(row+1, col) + " " + find(row, col+1));
        if (!find(row, col+1))
            return ADVANCE_THIS;
        if (!find(row+1, col))
            return ADVANCE_OTHER;
    }
    return ((bestRow==row)? ADVANCE_THIS: 0) |
        ((bestCol==col)? ADVANCE_OTHER: 0);

} // getExpandDirection()
	
unsigned char
Finder::getDistance(int row, int col) 
{
    if (find(row, col)) {
        return pm1->distance[row][col - pm1->first[row]];
    }
    std::cerr << "getDistance(" << row << "," << col << "): out of bounds" << std::endl;
    throw "getDistance index out of bounds";
} // getDistance()/2

void
Finder::setDistance(int row, int col, unsigned char b)
{
    if (find(row, col)) {
        pm1->distance[row][col - pm1->first[row]] = b;
        return;
    }
    std::cerr << "setDistance(" << row << "," << col << "," << b << "): out of bounds" << std::endl;
    throw "setDistance index out of bounds";
} // setDistance()

int
Finder::getPathCost(int row, int col)
{
    if (find(row, col)) // "1" avoids div by 0 below
        return pm1->bestPathCost[row][col - pm1->first[row]]*100/ (1+row+col);
    std::cerr << "getPathCost(" << row << "," << col << "): out of bounds" << std::endl;
    throw "getPathCost index out of bounds";
} // getPathCost()
	
int
Finder::getRawPathCost(int row, int col)
{
    if (find(row, col))
        return pm1->bestPathCost[row][col - pm1->first[row]];
    std::cerr << "getRawPathCost(" << row << "," << col << "): out of bounds" << std::endl;
    throw "getRawPathCost index out of bounds";
} // getRawPathCost()

void
Finder::setPathCost(int row, int col, int i)
{
    if (find(row, col)) {
         pm1->bestPathCost[row][col - pm1->first[row]] = i;
         return;
    }
    std::cerr << "setPathCost(" << row << "," << col << "," << i << "): out of bounds" << std::endl;
    throw "setPathCost index out of bounds";
} // setPathCost()

unsigned char
Finder::getDistance() 
{
    return pm1->distance[index1][index2];
} // getDistance()/0

void
Finder::setDistance(int b)
{
    pm1->distance[index1][index2] = (unsigned char)b;
} // setDistance()

int
Finder::getPathCost()
{
    return pm1->bestPathCost[index1][index2];
} // getPathCost()

void
Finder::setPathCost(int i)
{
    pm1->bestPathCost[index1][index2] = i;
} // setPathCost()

void
Finder::recalculatePathCostMatrix(int r1, int c1, int r2, int c2) 
{
    if (!find(r1,c1)) {
        std::cerr << "recalculatePathCostMatrix(" << r1 << "," << c1 << "): out of bounds" << std::endl;
        throw "recalculatePathCostMatrix index out of bounds";
    }
    int thisRowStart, c;
    int prevRowStart = 0, prevRowStop = 0;
    for (int r = r1; r <= r2; r++) {
        thisRowStart = pm1->first[r];
        if (thisRowStart < c1)
            thisRowStart = c1;
        for (c = thisRowStart; c <= c2; c++) {
            if (find(r,c)) {
                int i2 = index2;
                int newCost = pm1->distance[r][i2];
                int dir = 0;
                if (r > r1) {	// not first row
                    int min = -1;
                    if ((c > prevRowStart) && (c <= prevRowStop)) {
                        // diagonal from (r-1,c-1)
                        min = pm1->bestPathCost[r-1][c-pm1->first[r-1]-1] +
                            newCost * 2;
                        dir = ADVANCE_BOTH;
                    }
                    if ((c >= prevRowStart) && (c < prevRowStop)) {
                        // vertical from (r-1,c)
                        int cost = pm1->bestPathCost[r-1][c-pm1->first[r-1]] +
                            newCost;
                        if ((min == -1) || (cost < min)) {
                            min = cost;
                            dir = ADVANCE_THIS;
                        }
                    }
                    if (c > thisRowStart) {
                        // horizontal from (r,c-1)
                        int cost =pm1->bestPathCost[r][i2-1]+newCost;
                        if ((min == -1) || (cost < min)) {
                            min = cost;
                            dir = ADVANCE_OTHER;
                        }
                    }
                    pm1->bestPathCost[r][i2] = min;
                } else if (c > thisRowStart) {	// first row
                    // horizontal from (r,c-1)
                    pm1->bestPathCost[r][i2] = pm1->bestPathCost[r][i2-1] +
                        newCost;
                    dir = ADVANCE_OTHER;
                }
                if ((r != r1) || (c != c1)) {
                    pm1->distance[r][i2] = (unsigned char)
                        ((pm1->distance[r][i2] & MASK) | dir);
                }
            } else
                break;	// end of row
        }
        prevRowStart = thisRowStart;
        prevRowStop = c;
    }
} // recalculatePathCostMatrix()

int
Finder::retrievePath(bool smooth, vector<int> &pathx, vector<int> &pathy)
{
    int x = pm2->getFrameCount() - 1;
    int y = pm1->getFrameCount() - 1;

    pathx.clear();
    pathy.clear();

    while (find(y, x) && ((x > 0) || (y > 0))) {

        pathx.push_back(x);
        pathy.push_back(y);

        switch (getDistance() & ADVANCE_BOTH) {
        case ADVANCE_THIS:  y--; break;
        case ADVANCE_OTHER: x--; break;
        case ADVANCE_BOTH:  x--; y--; break;
        default: // this would indicate a bug, but we wouldn't want to hang
            cerr << "WARNING: Neither matcher advanced in path backtrack at (" << x << "," << y << ")" << endl;
            if (x > y) x--; else y--; break;
        }
    }

    std::reverse(pathx.begin(), pathx.end());
    std::reverse(pathy.begin(), pathy.end());

    if (smooth) {
        int smoothedLen = Path().smooth(pathx, pathy, pathx.size());
        return smoothedLen;
    } else {
        return pathx.size();
    }
}


