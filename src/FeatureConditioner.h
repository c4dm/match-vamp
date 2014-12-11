/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp feature extraction plugin using the MATCH audio alignment
    algorithm.

    Centre for Digital Music, Queen Mary, University of London.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef FEATURE_CONDITIONER_H
#define FEATURE_CONDITIONER_H

#include <vector>

/**
 * Take a series of feature vectors and apply conditioning of some
 * sort, such as normalisation or first-order derivative.
 *
 * Note that FeatureConditioner maintains internal frame-to-frame
 * state: use one FeatureConditioner per audio source, and construct a
 * new one for each new source.
 */
class FeatureConditioner
{
public:
    enum Normalisation {

        /** Do not normalise */
        NoNormalisation,
        
        /** Normalise each feature vector to have a sum of 1 */
        NormaliseToSum1,
        
        /** Normalise each feature vector by the long-term average of
         *  the summed energy */
        NormaliseToLTAverage,
    };

    enum OutputOrder {

	/** Output the normalised features without further processing */
	OutputFeatures,

	/** Output the half-wave rectified difference between the
	 *  previous and current features instead of the straight
	 *  feature values. */
	OutputRectifiedDerivative,

	/** Output the absolute difference between the previous and
	 *  current features instead of the straight feature
	 *  values. */
	OutputDerivative,
    };

    struct Parameters {

	Parameters() :
	    norm(NormaliseToSum1),
	    order(OutputRectifiedDerivative),
	    silenceThreshold(0.01),
	    decay(0.99)
	{}

	/** Feature normalisation. */
	Normalisation norm;

	/** Type of output to generate (plain feature, derivative etc). */
	OutputOrder order;

	/** Silence threshold. If non-zero, any feature whose total
	 *  energy (simply the sum of feature values) is below that
	 *  threshold will be rounded down to all zeros. Note that
	 *  this refers to the energy of the feature, not its
	 *  derivative even if that is what is being returned. */
	double silenceThreshold;

        /** Frame-to-frame decay factor in calculating long-term average */
	double decay;
    };
	
    /**
     * Construct a FeatureConditioner with the given parameters.
     *
     * Note that FeatureConditioner maintains internal frame-to-frame
     * state: use one FeatureConditioner per audio source, and construct
     * a new one for each new source.
     */
    FeatureConditioner(Parameters parameters) : m_params(parameters) { }

    /**
     * Process the given feature and return the conditioned feature.
     */
    std::vector<double> process(const std::vector<double> &feature);

protected:
    Parameters m_params;
    
    /** Long term average feature energy. */
    double m_ltAverage;

    /** The most recent feature, used for calculating the feature to
     *  feature difference. This is therefore not yet normalised. */
    std::vector<double> m_prev;
};

#endif
