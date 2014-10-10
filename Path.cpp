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

#include "Path.h"

int
Path::smooth(std::vector<int> &x, std::vector<int> &y, int length)
{
    if (length == 0)
        return 0;
    while ((int)val.size() < length) {
        val.push_back(0);
        len.push_back(0);
    }
    int p = 0;
    val[0] = len[0] = 0;
    for (int i = 1; i < length; i++) {	// H = 1; V = 2; D = 3
        int current = x[i] - x[i-1] + 2 * (y[i] - y[i-1]);
        if (current == val[p]) {
            len[p]++;
        } else if ((current == 3) || (val[p] == 0)) {
            val[++p] = current;
            len[p] = 1;
        } else if (val[p] + current == 3) {	// 1 + 2
            if (--len[p] == 0)
                p--;
            if (val[p] == 3)
                len[p]++;
            else {
                val[++p] = 3;
                len[p] = 1;
            }
        } else {	// val[p] == 3 && current != 3
            if ((val[p-1] == current) ||
                (val[p-1] == 0) ||
                (len[p] > MAX_RUN_LENGTH)) {
                val[++p] = current;
                len[p] = 1;
            } else {
                if (--len[p-1] == 0) {
                    val[p-1] = val[p];
                    len[p-1] = len[p];
                    p--;
                    if (val[p-1] == 3) {
                        len[p-1] += len[p];
                        p--;
                    }
                }
                len[p]++;
            }
        }
    }
    int i = 1;
    for (int pp = 1; pp <= p; pp++) {
        int dx = val[pp] & 1;
        int dy = val[pp] >> 1;
        for (int j = len[pp]; j > 0; j--, i++) {
            x[i] = x[i-1] + dx;
            y[i] = y[i-1] + dy;
        }
    }
    return i;
}
