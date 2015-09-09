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

#ifndef _PATH_H_
#define _PATH_H_

#include <vector>

class Path
{
public:
    Path() { }

    /** Smooths an alignment path.<BR>
     *  Consider the path as a sequence of horizontal (H), vertical (V) and
     *  diagonal (D) steps.  The smoothing consists of 2 rewrite rules:<BR>
     *  HnDmVn / Dm+n (where m is less than MAX_RUN_LENGTH)<BR>
     *  VnDmHn / Dm+n (where m is less than MAX_RUN_LENGTH)<BR>
     *  The new path is written over the old path.  Note that the end points of
     *  each application of a rewrite rule do not change.
     *  @return the length of the new path
     */
    int smooth(std::vector<int> &x, std::vector<int> &y, int length);

protected:
    static const int MAX_RUN_LENGTH = 50;

    std::vector<int> m_val;
    std::vector<int> m_len;
};

#endif

