/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "DistanceMetric.h"

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

static feature_t getTestFeature(double m)
{
    feature_t f;
    int fd[] = { 0, 1, 2, 3 };
    for (int i = 0; i < 4; ++i) {
        f.push_back(featurebin_t(fd[i] * m));
    }
    return f;
}

BOOST_AUTO_TEST_SUITE(TestDistanceMetric)

BOOST_AUTO_TEST_CASE(scale)
{
    DistanceMetric::Parameters params;
    params.scale = 1.0;
    DistanceMetric dm(params);
    
    BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(0.0), 0);
    BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(1.0), 1);
    BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(2.0), 2);

    if (sizeof(distance_t) == 1) {
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(256.0), 255);
    } else {
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(256.0), 256);
    }

    params.scale = 2.0;
    dm = DistanceMetric(params);
    
    if (sizeof(distance_t) == 1) {
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(0.0), 0);
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(1.0), 2);
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(2.0), 4);
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(128.0), 255);
    } else {
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(0.0), 0);
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(1.0), 1);
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(2.0), 2);
        BOOST_CHECK_EQUAL(dm.scaleValueIntoDistanceRange(256.0), 256);
    }
}

BOOST_AUTO_TEST_CASE(nonorm)
{
    DistanceMetric::Parameters params;
    params.norm = DistanceMetric::NoDistanceNormalisation;
    DistanceMetric dm(params);
    feature_t
        e1 = getTestFeature(1),
        e2 = getTestFeature(2),
        e0 = getTestFeature(0);

    double noise = 1e-3 * 4;
    
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), dm.scaleValueIntoDistanceRange(0.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), dm.scaleValueIntoDistanceRange(6.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), dm.scaleValueIntoDistanceRange(6.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), dm.scaleValueIntoDistanceRange(0.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), dm.scaleValueIntoDistanceRange(6.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), dm.scaleValueIntoDistanceRange(12.0 + noise));
}

BOOST_AUTO_TEST_CASE(sum)
{
    DistanceMetric::Parameters params;
    params.norm = DistanceMetric::NormaliseDistanceToSum;
    DistanceMetric dm(params);
    feature_t
        e1 = getTestFeature(1),
        e2 = getTestFeature(2),
        e0 = getTestFeature(0);

    double noise = 1e-3 * 4;
    
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), dm.scaleValueIntoDistanceRange(1.0));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), dm.scaleValueIntoDistanceRange(1.0));
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), dm.scaleValueIntoDistanceRange(1.0));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), dm.scaleValueIntoDistanceRange(noise / (12.0 + noise)));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), dm.scaleValueIntoDistanceRange((6.0 + noise) / (18.0 + noise)));
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), dm.scaleValueIntoDistanceRange(1.0));
}

BOOST_AUTO_TEST_SUITE_END()


    
