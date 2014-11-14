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
    if (!p1->m_firstPM)
        std::cerr << "Warning: wrong args in Finder()" << std::endl;
    pm1 = p1;
    pm2 = p2;
    index1 = 0;
    index2 = 0;
    rowRange = new int[2];
    colRange = new int[2];
    duration1 = -1;
    duration2 = -1;
} // constructor

Finder::~Finder()
{
    delete[] rowRange;
    delete[] colRange;
}

void
Finder::setDurations(int d1, int d2)
{
    duration1 = d1;
    duration2 = d2;
}

bool
Finder::find(int i1, int i2)
{
    if (i1 >= 0) {
        index1 = i1;
        index2 = i2 - pm1->m_first[i1];
    }
    return (i1 >= 0) && (i2 >= pm1->m_first[i1]) && (i2 < pm1->m_last[i1]);
} // find()

void
Finder::getColRange(int row, int *range)
{
    range[0] = pm1->m_first[row];
    range[1] = pm1->m_last[row];
} // getColRange()

void
Finder::getRowRange(int col, int *range)
{
    range[0] = pm2->m_first[col];
    range[1] = pm2->m_last[col];
} // getRowRange()

Matcher::Advance
Finder::getExpandDirection(int row, int col)
{
    return getExpandDirection(row, col, false);
} // getExpandDirection()

Matcher::Advance
Finder::getExpandDirection(int row, int col, bool check)
{
    double min = getPathCost(row, col);
    bestRow = row;
    bestCol = col;
    getRowRange(col, rowRange);
    if (rowRange[1] > row+1)
        rowRange[1] = row+1;	// don't cheat by looking at future :)
    for (int index = rowRange[0]; index < rowRange[1]; index++) {
        double tmp = getPathCost(index, col);
        if (tmp < min) {
            min = tmp;
            bestRow = index;
        }
    }
    getColRange(row, colRange);
    if (colRange[1] > col+1)
        colRange[1] = col+1;	// don't cheat by looking at future :)
    for (int index = colRange[0]; index < colRange[1]; index++) {
        double tmp = getPathCost(row, index);
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
        if (!find(row, col+1)) {
            return Matcher::AdvanceThis;
        } else if (!find(row+1, col)) {
            return Matcher::AdvanceOther;
        }
    }

    if (bestRow == row) {
        if (bestCol == col) {
            return Matcher::AdvanceBoth;
        } else {
            return Matcher::AdvanceThis;
        }
    } else if (bestCol == col) {
        return Matcher::AdvanceOther;
    } else {
        return Matcher::AdvanceNone;
    }

} // getExpandDirection()
	
float
Finder::getDistance(int row, int col) 
{
    if (find(row, col)) {
        return pm1->m_distance[row][col - pm1->m_first[row]];
    }
    std::cerr << "getDistance(" << row << "," << col << "): out of bounds" << std::endl;
    throw "getDistance index out of bounds";
} // getDistance()/2

void
Finder::setDistance(int row, int col, float b)
{
    if (find(row, col)) {
        pm1->m_distance[row][col - pm1->m_first[row]] = b;
        return;
    }
    std::cerr << "setDistance(" << row << "," << col << "," << b << "): out of bounds" << std::endl;
    throw "setDistance index out of bounds";
} // setDistance()

double
Finder::getPathCost(int row, int col)
{
    if (find(row, col)) // "1" avoids div by 0 below
        return pm1->m_bestPathCost[row][col - pm1->m_first[row]] / float(1+row+col);
    std::cerr << "getPathCost(" << row << "," << col << "): out of bounds" << std::endl;
    throw "getPathCost index out of bounds";
} // getPathCost()
	
double
Finder::getRawPathCost(int row, int col)
{
    if (find(row, col))
        return pm1->m_bestPathCost[row][col - pm1->m_first[row]];
    std::cerr << "getRawPathCost(" << row << "," << col << "): out of bounds" << std::endl;
    throw "getRawPathCost index out of bounds";
} // getRawPathCost()

void
Finder::setPathCost(int row, int col, double cost)
{
    if (find(row, col)) {
         pm1->m_bestPathCost[row][col - pm1->m_first[row]] = cost;
         return;
    }
    std::cerr << "setPathCost(" << row << "," << col << "," << cost << "): out of bounds" << std::endl;
    throw "setPathCost index out of bounds";
} // setPathCost()

Matcher::Advance
Finder::getAdvance()
{
    return pm1->m_advance[index1][index2];
}

void
Finder::setAdvance(Matcher::Advance a)
{
    pm1->m_advance[index1][index2] = a;
}

float
Finder::getDistance() 
{
    return pm1->m_distance[index1][index2];
} // getDistance()/0

void
Finder::setDistance(float b)
{
    pm1->m_distance[index1][index2] = b;
} // setDistance()

double
Finder::getPathCost()
{
    return pm1->m_bestPathCost[index1][index2];
} // getPathCost()

void
Finder::setPathCost(double cost)
{
    pm1->m_bestPathCost[index1][index2] = cost;
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
        thisRowStart = pm1->m_first[r];
        if (thisRowStart < c1)
            thisRowStart = c1;
        for (c = thisRowStart; c <= c2; c++) {
            if (find(r,c)) {
                int i2 = index2;
                float newCost = pm1->m_distance[r][i2];
                Matcher::Advance dir = Matcher::AdvanceNone;
                if (r > r1) {	// not first row
                    double min = -1;
                    if ((c > prevRowStart) && (c <= prevRowStop)) {
                        // diagonal from (r-1,c-1)
                        min = pm1->m_bestPathCost[r-1][c-pm1->m_first[r-1]-1] +
                            newCost * 2;
                        dir = Matcher::AdvanceBoth;
                    }
                    if ((c >= prevRowStart) && (c < prevRowStop)) {
                        // vertical from (r-1,c)
                        double cost = pm1->m_bestPathCost[r-1][c-pm1->m_first[r-1]] +
                            newCost;
                        if ((min == -1) || (cost < min)) {
                            min = cost;
                            dir = Matcher::AdvanceThis;
                        }
                    }
                    if (c > thisRowStart) {
                        // horizontal from (r,c-1)
                        double cost =pm1->m_bestPathCost[r][i2-1]+newCost;
                        if ((min == -1) || (cost < min)) {
                            min = cost;
                            dir = Matcher::AdvanceOther;
                        }
                    }
                    pm1->m_bestPathCost[r][i2] = min;
                } else if (c > thisRowStart) {	// first row
                    // horizontal from (r,c-1)
                    pm1->m_bestPathCost[r][i2] = pm1->m_bestPathCost[r][i2-1] +
                        newCost;
                    dir = Matcher::AdvanceOther;
                }
                pm1->m_advance[r][i2] = dir;
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

    if (duration2 > 0 && duration2 < pm2->getFrameCount()) {
        x = duration2 - 1;
    }
    if (duration1 > 0 && duration1 < pm1->getFrameCount()) {
        y = duration1 - 1;
    }

    recalculatePathCostMatrix(0, 0, y, x);

    pathx.clear();
    pathy.clear();

    while (find(y, x) && ((x > 0) || (y > 0))) {

//        cerr << "x = " << x << ", y = " << y;
        
        pathx.push_back(x);
        pathy.push_back(y);

        switch (getAdvance()) {
        case Matcher::AdvanceThis:
//            cerr << ", going down (dist = " << (int)getDistance() << ")" << endl;
            y--;
            break;
        case Matcher::AdvanceOther:
//            cerr << ", going left (dist = " << (int)getDistance() << ")" << endl;
            x--;
            break;
        case Matcher::AdvanceBoth:
//            cerr << ", going diag (dist = " << (int)getDistance() << ")" << endl;
            x--;
            y--;
            break;
        case Matcher::AdvanceNone: // this would indicate a bug, but we wouldn't want to hang
//            cerr << "WARNING: Neither matcher advanced in path backtrack at (" << x << "," << y << ")" << endl;
            if (x > y) {
                x--;
            } else {
                y--;
            }
            break;
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


