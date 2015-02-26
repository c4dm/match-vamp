
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
    
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), 0.0 + noise);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), 6.0 + noise);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), 6.0 + noise);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), 0.0 + noise);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), 6.0 + noise);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), 12.0 + noise);
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
    
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), 1.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), 1.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), 1.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), noise / (12.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), (6.0 + noise) / (18.0 + noise));
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), 1.0);
}

BOOST_AUTO_TEST_SUITE_END()


    
