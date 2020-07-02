/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright (c) 2007-2020 Simon Dixon, Chris Cannam, and Queen Mary
    University of London, Copyright (c) 2014-2015 Tido GmbH.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "FullDTW.h"

#include <stdexcept>

FullDTW::FullDTW(Parameters params, DistanceMetric::Parameters dparams) :
    m_params(params),
    m_metric(dparams)
{
}

pathcost_t
FullDTW::choose(CostOption x, CostOption y, CostOption d) {
    if (x.present && y.present) {
        if (!d.present) {
            throw std::logic_error("if x & y both exist, so must diagonal");
        }
        return std::min(std::min(x.cost, y.cost), d.cost);
    } else if (x.present) {
        return x.cost;
    } else if (y.present) {
        return y.cost;
    } else {
        return 0.0;
    }
}

pathcostmat_t
FullDTW::costSequences(const featureseq_t &s1, const featureseq_t &s2) {

    pathcostmat_t costs(s1.size(), pathcostvec_t(s2.size(), 0));

    for (size_t j = 0; j < s1.size(); ++j) {
        for (size_t i = 0; i < s2.size(); ++i) {
            distance_t dist = m_metric.calcDistance(s1[j], s2[i]);
            if (i == 0 && m_params.subsequence) {
                costs[j][i] = dist;
            } else {
                costs[j][i] = choose
                    (
                        { j > 0,
                          j > 0 ? dist + costs[j-1][i] : 0
                        },
                        { i > 0,
                          i > 0 ? dist + costs[j][i-1] : 0
                        },
                        { j > 0 && i > 0,
                          j > 0 && i > 0 ? dist + costs[j-1][i-1] : 0
                        });
            }
        }
    }

    return costs;
}

std::vector<size_t>
FullDTW::align(const featureseq_t &s1, const featureseq_t &s2) {

    // Return the index into s1 for each element in s2
        
    std::vector<size_t> alignment(s2.size(), 0);

    if (s1.empty() || s2.empty()) {
        return alignment;
    }

    auto costs = costSequences(s1, s2);

#ifdef DEBUG_DTW
    SVCERR << "Cost matrix:" << endl;
    for (auto v: cost) {
        for (auto x: v) {
            SVCERR << x << " ";
        }
        SVCERR << "\n";
    }
#endif

    size_t j = s1.size() - 1;
    size_t i = s2.size() - 1;

    if (m_params.subsequence) {
        pathcost_t min = 0.0;
        size_t minidx = 0;
        for (size_t j = 0; j < s1.size(); ++j) {
            if (j == 0 || costs[j][i] < min) {
                min = costs[j][i];
                minidx = j;
            }
        }
        j = minidx;
#ifdef DEBUG_DTW
        SVCERR << "Lowest cost at end of subsequence = " << min
               << " at index " << j << ", tracking back from there" << endl;
#endif
    }
        
    while (i > 0 && (j > 0 || m_params.subsequence)) {

        alignment[i] = j;
            
        pathcost_t a = costs[j-1][i];
        pathcost_t b = costs[j][i-1];
        pathcost_t both = costs[j-1][i-1];

        if (a < b) {
            --j;
            if (both <= a) {
                --i;
            }
        } else {
            --i;
            if (both <= b) {
                --j;
            }
        }
    }

    if (m_params.subsequence) {
        alignment[0] = j;
    }
        
    return alignment;
}
