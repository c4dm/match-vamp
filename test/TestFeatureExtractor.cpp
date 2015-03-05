
#include "FeatureExtractor.h"

#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

static int freq2mid(double freq)
{
    return int(round(57.0 + 12.0 * log(freq / 220.) / log(2.)));
}

static int freq2chroma(double freq)
{
    return freq2mid(freq) % 12;
}

static int bin2warped(int bin, int rate, int sz)
{
    // see comments in nonchroma below
    if (bin <= 33) return bin;
    double freq = (double(bin) * rate) / sz;
    int mid = freq2mid(freq);
    if (mid > 127) mid = 127;
    int outbin = mid - 77 + 33;
    return outbin;
}

BOOST_AUTO_TEST_SUITE(TestFeatureExtractor)

void checkAlternateProcessType(FeatureExtractor &fe, feature_t expected,
			       vector<double> real, vector<double> imag)
{
    vector<float> in;
    for (size_t i = 0; i < real.size(); ++i) {
	in.push_back(float(real[i]));
	in.push_back(float(imag[i]));
    }
    feature_t alt = fe.process(&in[0]);
    BOOST_CHECK_EQUAL_COLLECTIONS(alt.begin(), alt.end(),
				  expected.begin(), expected.end());
}

BOOST_AUTO_TEST_CASE(chroma)
{
    int szs[] = { 1024, 2048, 4000 };
    int rates[] = { 44100, 48000 };
    
    for (int irate = 0; irate < int(sizeof(rates)/sizeof(rates[0])); ++irate) {

	int rate = rates[irate];

	for (int isz = 0; isz < int(sizeof(szs)/sizeof(szs[0])); ++isz) {
	    
	    int sz = szs[isz];

	    int hs = sz / 2 + 1;
	    int fsz = 13;
    
	    FeatureExtractor::Parameters params { float(rate) };
	    params.fftSize = sz;
	    params.useChromaFrequencyMap = true;
	    params.minFrequency = 0;
	    FeatureExtractor fe(params);
	    BOOST_CHECK_EQUAL(fe.getFeatureSize(), fsz);

	    for (int bin = 0; bin < hs; ++bin) {

		vector<double> real, imag;
		real.resize(hs, 0.0);
		imag.resize(hs, 0.0);
	
		real[bin] += 10.0;
		imag[bin] += 10.0;

		// use two input sweeps, so we can test that they are
		// properly summed into the output bin
		real[hs-bin-1] += 5.0;
		imag[hs-bin-1] += 5.0;
	
		feature_t out = fe.process(real, imag);

		checkAlternateProcessType(fe, out, real, imag);

		// We expect to find all bins are 0 except for:
		// 
		// * two bins of 200 and 50 respectively, if the two input
		// bins are distinct and their output chroma are also distinct
		//
		// * one bin of value 250 (= 10^2 + 5^2), if the two input
		// bins are distinct but their output chroma are not
		//
		// * one bin of value 450 (= 15^2 + 15^2), if the input bins
		// are not distinct (the feature extractor sums energies
		// rather than magnitudes so as to integrate energy for a
		// partial in the face of spectral leakage)
		// 
		// The bin corresponding to each input frequency is
		// that of its semitone value (with C=0), except that
		// input bins less than the 17th are shepherded into
		// the separate bin 0 (see docs in FeatureExtractor.h)

		double cutoff = (17.0 * rate) / sz;
		
		feature_t expected(fsz);

		double infreq1 = (double(bin) * rate) / sz;

		if (bin == hs-bin-1) {
		    expected[freq2chroma(infreq1) + 1] += 450;
		} else {
		    if (infreq1 < cutoff) {
			expected[0] += 200;
		    } else {
			expected[freq2chroma(infreq1) + 1] += 200;
		    }
		    double infreq2 = (double(hs-bin-1) * rate) / sz;
		    if (infreq2 < cutoff) {
			expected[0] += 50;
		    } else {
			expected[freq2chroma(infreq2) + 1] += 50;
		    }
		}

		BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(),
					      expected.begin(), expected.end());
	    }
	}
    }
}

BOOST_AUTO_TEST_CASE(nonchroma)
{
    int rate = 44100;
    int sz = 2048;
    int hs = sz / 2 + 1;
    int fsz = 84;
    
    FeatureExtractor::Parameters params { float(rate) };
    params.fftSize = sz;
    params.minFrequency = 0;
    FeatureExtractor fe(params);
    BOOST_CHECK_EQUAL(fe.getFeatureSize(), fsz);

    for (int bin = 0; bin < hs; ++bin) {

	vector<double> real, imag;
	real.resize(hs, 0.0);
	imag.resize(hs, 0.0);
	
	real[bin] += 10.0;
	imag[bin] += 10.0;

	// use two input sweeps, so we can test that they are properly
	// summed into the output bin
	real[hs-bin-1] += 5.0;
	imag[hs-bin-1] += 5.0;

	feature_t out = fe.process(real, imag);

	checkAlternateProcessType(fe, out, real, imag);

	// We expect to find all bins are 0 except for:
	// 
	// * two bins of 200 and 50 respectively, if the two input
	// bins are distinct and their output bins are also distinct
	//
	// * one bin of value 250 (= 10^2 + 5^2), if the two input
	// bins are distinct but their output bins are not
	//
	// * one bin of value 450 (= 15^2 + 15^2), if the input bins
	// are not distinct (the feature extractor sums energies
	// rather than magnitudes so as to integrate energy for a
	// partial in the face of spectral leakage)
	//
	// The first 34 input bins (i.e. up to and including bin 33,
	// 733Hz, MIDI pitch 77.something) are mapped linearly, those
	// above that and up to MIDI pitch 127 (12544Hz) are mapped
	// logarithmically, remaining bins are all mapped into the
	// final output bin.
	//
	// So MIDI pitches up to and including 77 are mapped linearly
	// by frequency into 34 bins; those from 78-126 inclusive are
	// mapped linearly by MIDI pitch into the next 49 bins;
	// everything above goes into the last bin, for 84 bins total.

	feature_t expected(fsz);

	if (bin == hs-bin-1) {
	    expected[bin2warped(bin, rate, sz)] += 450;
	} else {
	    expected[bin2warped(bin, rate, sz)] += 200;
	    expected[bin2warped(hs-bin-1, rate, sz)] += 50;
	}

	BOOST_CHECK_EQUAL_COLLECTIONS(out.begin(), out.end(),
				      expected.begin(), expected.end());
    }
}

BOOST_AUTO_TEST_SUITE_END()
