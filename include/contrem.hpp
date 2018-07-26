/*
 * processor.hpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#ifndef PROCESSOR_HPP_
#define PROCESSOR_HPP_

#include "reader.hpp"

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
	std::string outfile;	///<! The output file.
	std::string driver;		///<! The output driver.
	std::string extension;	///<! The output extension.
	int bufferSize;			///<! The buffer size.
	int threads;			///<! The number of threads to use.
	bool sampleStats;		///<! True if sample stats should be used, otherwise population stats.

	/**
	 * Process the continuum removal job.
	 *
	 * @param listener An implementation of ContremListener that can receive status events.
	 * @param reader A Reader implementation that can provide data.
	 * @param config A configuration object.
	 */
	void run(ContremListener* listener, Reader* reader);

	/**
	 * Returns the current processing progress as a double from 0 to 1.
	 *
	 * @return The current processing progress as a double from 0 to 1.
	 */
	double progress() const;

};

#endif /* PROCESSOR_HPP_ */
