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
#include <iomanip>

using namespace std;

Finder::Finder(Matcher *pm)
{
    m_m = pm;
    m_duration1 = -1;
    m_duration2 = -1;
} // constructor

Finder::~Finder()
{
}

void
Finder::setDurations(int d1, int d2)
{
    m_duration1 = d1;
    m_duration2 = d2;
}

Matcher::Advance
Finder::getExpandDirection(int row, int col)
{
    double min = m_m->getPathCost(row, col);
    
    int bestRow = row;
    int bestCol = col;

    pair<int, int> rowRange = m_m->getRowRange(col);
    if (rowRange.second > row+1) {
        rowRange.second = row+1;	// don't cheat by looking at future :)
    }
    for (int index = rowRange.first; index < rowRange.second; index++) {
        double tmp = m_m->getPathCost(index, col);
        if (tmp < min) {
            min = tmp;
            bestRow = index;
        }
    }

    pair<int, int> colRange = m_m->getColRange(row);
    if (colRange.second > col+1) {
        colRange.second = col+1;	// don't cheat by looking at future :)
    }
    for (int index = colRange.first; index < colRange.second; index++) {
        double tmp = m_m->getPathCost(row, index);
        if (tmp < min) {
            min = tmp;
            bestCol = index;
            bestRow = row;
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
}

void
Finder::recalculatePathCostMatrix(int r1, int c1, int r2, int c2)
{
    recalculatePathCostMatrix(r1, c1, r2, c2, m_m->getDiagonalWeight());
}

void
Finder::recalculatePathCostMatrix(int r1, int c1, int r2, int c2,
                                  float diagonalWeight) 
{
    int prevRowStart = 0, prevRowStop = 0;
    
    for (int r = r1; r <= r2; r++) {

        pair<int, int> colRange = m_m->getColRange(r);

        int rowStart = max(c1, colRange.first);
        int rowStop = min(c2 + 1, colRange.second);
        
        for (int c = rowStart; c < rowStop; c++) {

            float newCost = m_m->getDistance(r, c);
            Matcher::Advance dir = Matcher::AdvanceNone;

            if (r > r1) {	// not first row
                double min = -1;
                if ((c > prevRowStart) && (c <= prevRowStop)) {
                    // diagonal from (r-1,c-1)
                    min = m_m->getPathCost(r-1, c-1) + newCost * diagonalWeight;
                    dir = Matcher::AdvanceBoth;
                }
                if ((c >= prevRowStart) && (c < prevRowStop)) {
                    // vertical from (r-1,c)
                    double cost = m_m->getPathCost(r-1, c) + newCost;
                    if ((min < 0) || (cost < min)) {
                        min = cost;
                        dir = Matcher::AdvanceThis;
                    }
                }
                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    double cost = m_m->getPathCost(r, c-1) + newCost;
                    if ((min < 0) || (cost < min)) {
                        min = cost;
                        dir = Matcher::AdvanceOther;
                    }
                }
                
                m_m->setPathCost(r, c, dir, min);

            } else if (c > rowStart) {	// first row
                // horizontal from (r,c-1)
                m_m->setPathCost(r, c, Matcher::AdvanceOther,
                                   m_m->getPathCost(r, c-1) + newCost);
            }
        }

        prevRowStart = rowStart;
        prevRowStop = rowStop;
    }
} 

#ifdef PERFORM_ERROR_CHECKS
Finder::ErrorPosition
Finder::checkPathCostMatrix() 
{
    ErrorPosition err;

    int r1 = 0;
    int c1 = 0;
    int r2 = m_m->getFrameCount() - 1;
    int c2 = m_m->getOtherFrameCount() - 1;

    if (r2 < r1 || c2 < c1) {
        return err;
    }

    int prevRowStart = 0, prevRowStop = 0;

    float diagonalWeight = m_m->getDiagonalWeight();
    
    for (int r = r1; r <= r2; r++) {

        pair<int, int> colRange = m_m->getColRange(r);

        int rowStart = max(c1, colRange.first);
        int rowStop = min(c2 + 1, colRange.second);
        
        for (int c = rowStart; c < rowStop; c++) {

            float newCost = m_m->getDistance(r, c);
            double updateTo = -1.0;
            Matcher::Advance dir = Matcher::AdvanceNone;

            if (r > r1) { // not first row
                double min = -1;
                if ((c > prevRowStart) && (c <= prevRowStop)) {
                    // diagonal from (r-1,c-1)
                    min = m_m->getPathCost(r-1, c-1) + newCost * diagonalWeight;
                    err.prevCost = m_m->getPathCost(r-1, c-1);
                    err.distance = newCost * diagonalWeight;
                    dir = Matcher::AdvanceBoth;
                }
                if ((c >= prevRowStart) && (c < prevRowStop)) {
                    // vertical from (r-1,c)
                    double cost = m_m->getPathCost(r-1, c) + newCost;
                    if ((min < 0) || (cost < min)) {
                        min = cost;
                        err.prevCost = m_m->getPathCost(r-1, c);
                        err.distance = newCost;
                        dir = Matcher::AdvanceThis;
                    }
                }
                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    double cost = m_m->getPathCost(r, c-1) + newCost;
                    if ((min < 0) || (cost < min)) {
                        min = cost;
                        err.prevCost = m_m->getPathCost(r, c-1);
                        err.distance = newCost;
                        dir = Matcher::AdvanceOther;
                    }
                }

                updateTo = min;

            } else { // first row

                if (c > rowStart) {
                    // horizontal from (r,c-1)
                    updateTo = m_m->getPathCost(r, c-1) + newCost;
                    err.prevCost = m_m->getPathCost(r, c-1);
                    err.distance = newCost;
                    dir = Matcher::AdvanceOther;
                }
            }

            if (dir != Matcher::AdvanceNone) {
                if (m_m->getAdvance(r, c) != dir) {
                    err.type = ErrorPosition::WrongAdvance;
                    err.r = r;
                    err.c = c;
                    err.costWas = m_m->getPathCost(r, c);
                    err.costShouldBe = updateTo;
                    err.advanceWas = m_m->getAdvance(r, c);
                    err.advanceShouldBe = dir;
                    return err;
                }
                if (m_m->getPathCost(r, c) != updateTo) {
                    err.type = ErrorPosition::WrongCost;
                    err.r = r;
                    err.c = c;
                    err.costWas = m_m->getPathCost(r, c);
                    err.costShouldBe = updateTo;
                    err.advanceWas = m_m->getAdvance(r, c);
                    err.advanceShouldBe = dir;
                    return err;
                }
            } else {
                // AdvanceNone should occur only at r = r1, c = c1
                if (r != r1 || c != c1) {
                    err.type = ErrorPosition::NoAdvance;
                    err.r = r;
                    err.c = c;
                    err.costWas = m_m->getPathCost(r, c);
                    err.costShouldBe = updateTo;
                    err.advanceWas = m_m->getAdvance(r, c);
                    err.advanceShouldBe = dir;
                    return err;
                }
            }
        }

        prevRowStart = rowStart;
        prevRowStop = rowStop;
    }

    return err;
}

void
Finder::checkAndReport()
{
    cerr << "Finder: Checking path-cost matrix..." << endl;
    ErrorPosition err = checkPathCostMatrix();
    if (err.type == ErrorPosition::NoError) {
        cerr << "No errors found" << endl;
    } else {
        cerr << "\nWARNING: Checking path-cost matrix returned mismatch:" << endl;
        cerr << "Type: " << err.type << ": ";
        switch (err.type) {
        case ErrorPosition::NoError: break;
        case ErrorPosition::WrongCost: cerr << "WrongCost"; break;
        case ErrorPosition::WrongAdvance: cerr << "WrongAdvance"; break;
        case ErrorPosition::NoAdvance: cerr << "NoAdvance"; break;
        }
        cerr << endl;
        cerr << "At row " << err.r << ", column " << err.c
             << "\nShould be advancing "
             << Matcher::advanceToString(err.advanceShouldBe)
             << ", advance in matrix is "
             << Matcher::advanceToString(err.advanceWas)
             << "\nPrev cost " << err.prevCost
             << " plus distance " << err.distance << " gives "
             << err.costShouldBe << ", matrix contains " << err.costWas
             << endl;
        cerr << "Note: diagonal weight = " << m_m->getDiagonalWeight() << endl;
        cerr << endl;

        int w(4);
        int ww(15);

        cerr << "Distance matrix leading up to this point:" << endl;
        cerr << setprecision(12) << setw(w) << "";
        for (int i = -4; i <= 0; ++i) {
            cerr << setw(ww) << i;
        }
        cerr << endl;
        for (int j = -4; j <= 0; ++j) {
            cerr << setw(w) << j;
            for (int i = -4; i <= 0; ++i) {
                cerr << setw(ww) << m_m->getDistance(err.r + j, err.c + i);
            }
            cerr << endl;
        }
        cerr << endl;

        cerr << "Cost matrix leading up to this point:" << endl;
        cerr << setw(w) << "";
        for (int i = -4; i <= 0; ++i) {
            cerr << setw(ww) << i;
        }
        cerr << endl;
        for (int j = -4; j <= 0; ++j) {
            cerr << setw(w) << j;
            for (int i = -4; i <= 0; ++i) {
                cerr << setw(ww) << m_m->getPathCost(err.r + j, err.c + i);
            }
            cerr << endl;
        }
        cerr << endl;
    }
}
#endif

void
Finder::retrievePath(vector<int> &pathx,
                     vector<int> &pathy,
                     vector<float> &distances)
{
    pathx.clear();
    pathy.clear();
    distances.clear();

#ifdef PERFORM_ERROR_CHECKS
//    checkAndReport();
#endif

    int ex = m_m->getOtherFrameCount() - 1;
    int ey = m_m->getFrameCount() - 1;

    if (ex < 0 || ey < 0) {
        return;
    }
    
    int x = ex;
    int y = ey;

    if (m_duration2 > 0 && m_duration2 < m_m->getOtherFrameCount()) {
        x = m_duration2 - 1;
    }
    if (m_duration1 > 0 && m_duration1 < m_m->getFrameCount()) {
        y = m_duration1 - 1;
    }
    
//    cerr << "before: x = " << x << ", y = " << y << endl;

    if (!m_m->isAvailable(y, x)) {
        // Path did not pass through the expected end point --
        // probably means the pieces are substantially different in
        // the later bits. Reset the expected end point to the end of
        // both files including any trailing silence.
        cerr << "NOTE: Path did not pass through expected end point, inputs are probably significantly different" << endl;
        x = ex;
        y = ey;
    }

//    recalculatePathCostMatrix(0, 0, y, x);

//    cerr << "start: x = " << x << ", y = " << y << endl;
    
    while (m_m->isAvailable(y, x) && (x > 0 || y > 0)) {

//        cerr << "x = " << x << ", y = " << y;
        
        pathx.push_back(x);
        pathy.push_back(y);
        distances.push_back(m_m->getDistance(y, x));

        switch (m_m->getAdvance(y, x)) {
        case Matcher::AdvanceThis:
//            cerr << ", going down (dist = " << getDistance() << ")" << endl;
            y--;
            break;
        case Matcher::AdvanceOther:
//            cerr << ", going left (dist = " << getDistance() << ")" << endl;
            x--;
            break;
        case Matcher::AdvanceBoth:
//            cerr << ", going diag (dist = " << getDistance() << ")" << endl;
            x--;
            y--;
            break;
        case Matcher::AdvanceNone: // this would indicate a bug, but we wouldn't want to hang
            cerr << "WARNING: Neither matcher advanced in path backtrack at (" << x << "," << y << ")" << endl;
            if (x > y) {
                x--;
            } else {
                y--;
            }
            break;
        }
    }

    if (x > 0 || y > 0) {
        cerr << "WARNING: Ran out of available path at (" << y << "," << x
             << ")!" << endl;
    }
    
    reverse(pathx.begin(), pathx.end());
    reverse(pathy.begin(), pathy.end());
    reverse(distances.begin(), distances.end());
}

void
Finder::smoothWithPinPoints(const map<int, int> &pinpoints)
{
    cerr << "Pin points are:" << endl;

    typedef map<int, int> PPMap;

    for (PPMap::const_iterator i = pinpoints.begin(); i != pinpoints.end(); ++i) {
        cerr << "[" << i->first << "," << i->second << "] ";
    }
    
    cerr << endl;

    if (pinpoints.size() < 2) return;

    pair<int, int> prev = *pinpoints.begin();

    for (PPMap::const_iterator i = pinpoints.begin(); i != pinpoints.end(); ++i) {
        if (i == pinpoints.begin()) continue;

        pair<int, int> curr = *i;

        recalculatePathCostMatrix(prev.second, prev.first,
                                  curr.second, curr.first, 1.0);
        
        prev = curr;
    }
}


