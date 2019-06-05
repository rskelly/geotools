/*
 * processor.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef _CONTREM_HPP_
#define _CONTREM_HPP_

#include "reader.hpp"

namespace hlrg {

/**
 * Forward declaration.
 */
class Contrem;

/**
 * An interface that provides implementors with the ability to receive
 * notifications from the Contrem.
 */
class ContremListener {
public:

	/**
	 * Called when the convolution has started.
	 *
	 * @param A Contrem.
	 */
	virtual void started(Contrem*) = 0;

	/**
	 * Called when the convolution status has updated.
	 *
	 * @param A Contrem.
	 */
	virtual void update(Contrem*) = 0;

	/**
	 * Called when the convolution has stopped.
	 *
	 * @param A Contrem.
	 */
	virtual void stopped(Contrem*) = 0;

	/**
	 * Called when the convolution has finished.
	 *
	 * @param A Contrem.
	 */
	virtual void finished(Contrem*) = 0;

	virtual ~ContremListener() {}
};


/**
 * Performs the continuum removal process.
 */
class Contrem {
public:
	std::string output;			///<! The output file.
	std::string outputType;		///<! The output file type.
	std::string extension;		///<! The output extension.
	std::string roi;			///<! The mask/ROI; Shapefile, SQLite, ENVI ROI.
	std::string roiType;		///<! The mask/ROI file type. May be empty.
	std::string spectra;		///<! The input spectra; raster or CSV.
	std::string spectraType;	///<! The input spectra file type. May be empty.
	double minWl;				///<! The lower bound of the wavelength range to process.
	double maxWl;				///<! The upper bound of the wavelength range to process.
	int threads;				///<! The number of threads to use.
	bool sampleStats;			///<! True if sample stats should be used, otherwise population stats.

	/**
	 * Process the continuum removal job.
	 *
	 * @param listener An implementation of ContremListener that can receive status events.
	 * @param config A configuration object.
	 */
	void run(ContremListener* listener);

	/**
	 * Returns the current processing progress as a double from 0 to 1.
	 *
	 * @return The current processing progress as a double from 0 to 1.
	 */
	double progress() const;

};

} // hlrg

#endif /* _CONTREM_HPP_ */
