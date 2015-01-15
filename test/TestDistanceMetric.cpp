
#include "DistanceMetric.h"

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

static vector<double> getTestFeature(double m)
{
    vector<double> f;
    double fd[] = { 0, 1, 2, 3 };
    for (int i = 0; i < 4; ++i) {
	f.push_back(fd[i] * m);
    }
    return f;
}

BOOST_AUTO_TEST_SUITE(TestDistanceMetric)

BOOST_AUTO_TEST_CASE(nonorm)
{
    DistanceMetric::Parameters params;
    params.norm = DistanceMetric::NoDistanceNormalisation;
    params.silencePenalty = 0.0;
    DistanceMetric dm(params);
    vector<double>
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	e0 = getTestFeature(0);

    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), 0.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), 6.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), 6.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), 0.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), 6.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), 12.0);
}

BOOST_AUTO_TEST_CASE(sum)
{
    DistanceMetric::Parameters params;
    params.norm = DistanceMetric::NormaliseDistanceToSum;
    params.silencePenalty = 0.0;
    DistanceMetric dm(params);
    vector<double>
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	e0 = getTestFeature(0);

    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), 0.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), 1.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), 1.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), 0.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), 1.0/3.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), 1.0);
}

BOOST_AUTO_TEST_CASE(penalty)
{
    DistanceMetric::Parameters params;
    params.norm = DistanceMetric::NormaliseDistanceToSum;
    params.silencePenalty = 0.5;
    DistanceMetric dm(params);
    vector<double>
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	e0 = getTestFeature(0);

    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e0), 0.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e0), 1.5);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e1), 1.5);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e1), 0.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e1, e2), 1.0/3.0);
    BOOST_CHECK_EQUAL(dm.calcDistance(e0, e2), 1.5);
}

BOOST_AUTO_TEST_SUITE_END()


    
