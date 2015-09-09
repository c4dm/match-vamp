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

#include "Path.h"

int
Path::smooth(std::vector<int> &x, std::vector<int> &y, int length)
{
    if (length <= 0)
        return 0;
    while (m_val.size() < static_cast<std::vector<int>::size_type>(length)) {
        m_val.push_back(0);
        m_len.push_back(0);
    }
    int p = 0;
    m_val[0] = m_len[0] = 0;
    for (int i = 1; i < length; i++) {	// H = 1; V = 2; D = 3
        int current = x[i] - x[i-1] + 2 * (y[i] - y[i-1]);
        if (current == m_val[p]) {
            m_len[p]++;
        } else if ((current == 3) || (m_val[p] == 0)) {
            m_val[++p] = current;
            m_len[p] = 1;
        } else if (m_val[p] + current == 3) {	// 1 + 2
            if (--m_len[p] == 0)
                p--;
            if (m_val[p] == 3)
                m_len[p]++;
            else {
                m_val[++p] = 3;
                m_len[p] = 1;
            }
        } else {	// m_val[p] == 3 && current != 3
            if ((m_val[p-1] == current) ||
                (m_val[p-1] == 0) ||
                (m_len[p] > MAX_RUN_LENGTH)) {
                m_val[++p] = current;
                m_len[p] = 1;
            } else {
                if (--m_len[p-1] == 0) {
                    m_val[p-1] = m_val[p];
                    m_len[p-1] = m_len[p];
                    p--;
                    if (m_val[p-1] == 3) {
                        m_len[p-1] += m_len[p];
                        p--;
                    }
                }
                m_len[p]++;
            }
        }
    }
    int i = 1;
    for (int pp = 1; pp <= p; pp++) {
        int dx = m_val[pp] & 1;
        int dy = m_val[pp] >> 1;
        for (int j = m_len[pp]; j > 0; j--, i++) {
            x[i] = x[i-1] + dx;
            y[i] = y[i-1] + dy;
        }
    }
    return i;
}
