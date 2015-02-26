
#include "FeatureConditioner.h"

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

BOOST_AUTO_TEST_SUITE(TestFeatureConditioner)

BOOST_AUTO_TEST_CASE(nonorm_features)
{
    FeatureConditioner::Parameters params;
    params.norm = FeatureConditioner::NoNormalisation;
    params.order = FeatureConditioner::OutputFeatures;
    feature_t
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	e0 = getTestFeature(0);

    params.silenceThreshold = 1.0;
    FeatureConditioner fc(params);
    feature_t out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e2.begin(), e2.end());

    params.silenceThreshold = 7.0;
    fc = FeatureConditioner(params);
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e2.begin(), e2.end());
}
    
BOOST_AUTO_TEST_CASE(nonorm_rectderiv)
{
    FeatureConditioner::Parameters params;
    params.norm = FeatureConditioner::NoNormalisation;
    params.order = FeatureConditioner::OutputRectifiedDerivative;
    feature_t
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	e0 = getTestFeature(0);

    params.silenceThreshold = 1.0;
    FeatureConditioner fc(params);
    feature_t out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());

    params.silenceThreshold = 7.0;
    fc = FeatureConditioner(params);
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
}
    
BOOST_AUTO_TEST_CASE(nonorm_deriv)
{
    FeatureConditioner::Parameters params;
    params.norm = FeatureConditioner::NoNormalisation;
    params.order = FeatureConditioner::OutputDerivative;
    feature_t
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	e3 = getTestFeature(3),
	e0 = getTestFeature(0);

    params.silenceThreshold = 1.0;
    FeatureConditioner fc(params);
    feature_t out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());

    params.silenceThreshold = 7.0;
    fc = FeatureConditioner(params);
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e3);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e2.begin(), e2.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e1.begin(), e1.end());
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
}
    
BOOST_AUTO_TEST_CASE(sum1_features)
{
    FeatureConditioner::Parameters params;
    params.norm = FeatureConditioner::NormaliseToSum1;
    params.order = FeatureConditioner::OutputFeatures;
    feature_t
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	en = getTestFeature(1.0/6.0),
	e0 = getTestFeature(0);

    params.silenceThreshold = 1.0;
    FeatureConditioner fc(params);
    feature_t out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en.begin(), en.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en.begin(), en.end());

    params.silenceThreshold = 7.0;
    fc = FeatureConditioner(params);
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en.begin(), en.end());
}
    
BOOST_AUTO_TEST_CASE(sum1_rectderiv)
{
    FeatureConditioner::Parameters params;
    params.norm = FeatureConditioner::NormaliseToSum1;
    params.order = FeatureConditioner::OutputRectifiedDerivative;
    feature_t
	e1 = getTestFeature(1),
	e2 = getTestFeature(2),
	en = getTestFeature(1.0/6.0),
	en2 = getTestFeature(1.0/12.0),
	e0 = getTestFeature(0);

    params.silenceThreshold = 1.0;
    FeatureConditioner fc(params);
    feature_t out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en.begin(), en.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en2.begin(), en2.end());
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());

    params.silenceThreshold = 7.0;
    fc = FeatureConditioner(params);
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en2.begin(), en2.end());
    out = fc.process(e1);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), e0.begin(), e0.end());
    out = fc.process(e2);
    BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(), en2.begin(), en2.end());
}
    

BOOST_AUTO_TEST_SUITE_END()

