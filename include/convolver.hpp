/*
 * convolver.hpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#ifndef INCLUDE_CONVOLVER_HPP_
#define INCLUDE_CONVOLVER_HPP_

#include <string>
#include <fstream>

/**
 * The kernel class holds the step values, plus the min
 * and max wavelengths for the range.
 */
class Kernel {
private:
	double m_wl;
	double m_fwhm;
	double m_threshold;
	double m_sigma;
	double m_halfWidth;

public:

	/**
	 * Create a kernel around the given wavelength with the given
	 * full width at half maximum. The threshold gives the value
	 * of the Gaussian beyond which the kernel is not applied.
	 *
	 * @param wl The centre wavelength.
	 * @param fwhm The full width at half maximum.
	 * @param threshold The minimum function value.
	 */
	Kernel(double wl, double fwhm, double threshold);

	/**
	 * Set the wavelength.
	 *
	 * @param wl The wavelength.
	 */
	void setWavelength(double wl);

	/**
	 * Get the wavelength.
	 *
	 * @return The wavelength.
	 */
	double wl() const;

	/**
	 * Set the full width at half maximum.
	 *
	 * @param fwhm The full width at half maximum.
	 */
	void setFWHM(double fwhm);

	/**
	 * Get the full width at half maximum.
	 *
	 * @return The full width at half maximum.
	 */
	double fwhm() const;

	/**
	 * Initialize the internal variables.
	 */
	void init();

	/**
	 * Returns the absolute distance from the centre wavelength, to the wavelength where
	 * the function falls below the given threshold. This gives the wavelength
	 * range to which the kernel can be applied.
	 *
	 * @param threshold The minimum y-value of the function.
	 * @return The distance (in wavelength) from the centre of the
	 * 		   curve to the point where it falls below the threshold.
	 */
	double halfWidth() const;

	/**
	 * Return the value of the Gaussian for the given wavelength.
	 *
	 * @param The target wavelength.
	 * @return The Gaussian function output.
	 */
	double operator()(double wl0) const;

	/**
	 * Returns the threshold, the minimum value of the Gaussian
	 * before it is ignored.
	 *
	 * @return The minimum Gaussian value.
	 */
	double threshold() const;

};


/**
 * Contains the information about each band that will be used to build
 * the convolution kernel.
 */
class BandProp {
public:
	int band;		///<! The 1-based band index.
	double wl;		///<! The wavelength.
	double fwhm;	///<! Full width at half maximum.

	/**
	 * Create a BandProp instance.
	 *
	 * @param band The 1-based band index.
	 * @param wl The wavelength.
	 * @param fwhm The full width of half maximum.
	 */
	BandProp(int band, double wl, double fwhm);
};


/**
 * Represents a single band in the input/output spectra.
 */
class Band {
public:
	double wl;			///!< The wavelength.
	double value;		///!< The value or intensity.
	double scaledValue; ///!< The scaled value of the intensity. This is used for calculations. Is identical to value by default.

	/**
	 * Create a band.
	 *
	 * @param wl The wavelength.
	 * @param value The value or intensity.
	 */
	Band(double wl, double value);

	Band();
};



/**
 * The spectrum object has a list of bands. When loaded from a file,
 * loads the bands and stores the header properties.
 *
 * When the next() method is called, reads the next available spectrum,
 * updates the bands and the date and timestamp.
 */
class Spectrum {
private:
	std::ifstream m_input;							///<! The input stream.
	std::string m_buf;								///<! The current file buffer.
	size_t m_count;									///<! The number of records in the file.

public:
	std::vector<Band> bands;						///<! A list of the bands. This changes as the file is read through.
	std::map<std::string, std::string> properties;	///<! Properties read from the header block.
	std::string date;								///<! The date of the current row.
	long time;										///<! The timestamp of the current row.

	/**
	 * Load the data file and read the header information.
	 * After load is called (and returns true), the next method must be
	 * called to populate the bands list with data.
	 *
	 * @param filename A data file.
	 * @return True if the file is loaded and has information in it.
	 */
	bool load(const std::string& filename);

	/**
	 * Return the number of spectra.
	 */
	size_t count() const;

	/**
	 * Advance the reader to the next line of data. This becomes the Spectrum's current state:
	 * the list of Bands contains data from the current row.
	 *
	 * @return True if a row has been read, false if there were none left.
	 */
	bool next();

	/**
	 * Sets up the given Spectrum with the same structure and
	 * wavelengths as the current one.
	 *
	 * @param spec A Spectrume.
	 */
	void setup(Spectrum& spec);

	/**
	 * Apply this scaling factor to every band's intensity value.
	 *
	 * @param scale A scale factor.
	 */
	void scale(double scale);

	/**
	 * Convolve using the given Kernel, writing the output to the given
	 * Band.
	 *
	 * @param kernel A Kernel.
	 * @param band A Band.
	 */

	void convolve(Kernel& kernel, Band& band);

	/**
	 * Set the values of all bands to zero.
	 */
	void reset();

	/**
	 * Write the header of an output spectrum, between the given
	 * wavelengths.
	 *
	 * @param out An output stream.
	 * @param minWl The minimum wavelength to print.
	 * @param maxWl The maximum wavelength to print.
	 */
	void writeHeader(std::ostream& out, double minWl, double maxWl);

	/**
	 * Write the bands of an output spectrum, between the given
	 * wavelengths.
	 *
	 * @param out An output stream.
	 * @param minWl The minimum wavelength to print.
	 * @param maxWl The maximum wavelength to print.
	 */
	void write(std::ostream& out, double minWl, double maxWl);

};


/**
 * Read the convolution configuration file. The columns are
 * band #, wavelength, full width at half maximum.
 * Names are ignored; the order is important.
 */
class BandPropsReader {
public:
	std::map<int, BandProp> bandProps;	///<! The map of BandProp object.
	double minWl;						///<! The minimum wavelength in the band definition file.
	double maxWl;						///<! The maximum wavelength in the band definition file.

	/**
	 * Load the band definition file.
	 *
	 * @param filename The filename.
	 */
	void load(const std::string& filename);

	/**
	 * Returns a vector containing the list of band numbers.
	 *
	 * @return A vector containing the list of band numbers.
	 */
	std::vector<int> bands() const;

	/**
	 * Configure the given kernel using information for the given 1-based
	 * band index. Throws a runtime exception if the band is not available.
	 *
	 * @param kernel The Kernel.
	 * @param band The 1-based band index.
	 */
	void configureKernel(Kernel& kernel, int band);

	/**
	 * Configure the given Spectrum with the bands represented
	 * by the band definition file.
	 *
	 * @param spec A Spectrum object.
	 */
	void configureSpectrum(Spectrum& spec);

};



/**
 * Forward declaration.
 */
class Convolver;

/**
 * An interface that provides implementors with the ability to receive
 * notifications from the Convolver.
 */
class ConvolverListener {
public:

	/**
	 * Called when the convolution has started.
	 *
	 * @param A Convolver.
	 */
	virtual void started(Convolver*) = 0;

	/**
	 * Called when the convolution status has updated.
	 *
	 * @param A Convolver.
	 */
	virtual void update(Convolver*) = 0;

	/**
	 * Called when the convolution has stopped.
	 *
	 * @param A Convolver.
	 */
	virtual void stopped(Convolver*) = 0;

	/**
	 * Called when the convolution has finished.
	 *
	 * @param A Convolver.
	 */
	virtual void finished(Convolver*) = 0;

	virtual ~ConvolverListener() {}
};


/**
 * Performs the convolution of a spectral file using Gaussian kernels
 * parameterized from a band definition file.
 */
class Convolver {
private:
	double m_progress;	///<! The current progress, between 0 and 1.

public:

	/**
	 * Run the convolver on the given files. The listener will receive callbacks.
	 *
	 * @param listener A ConvolverListener to receive updates.
	 * @param bandDef The band definition file.
	 * @param spectra The spectral data file.
	 * @param output The output file.
	 * @param inputScale Scale every input spectral value by this much.
	 * @param running A reference to a boolean that is true so long as the processor should keep running.
	 */
	void run(ConvolverListener& listener,
			const std::string& bandDef, const std::string& spectra, const std::string& output,
			double inputScale, bool& running);

	/**
	 * Return the progress as a double between 0 and 1.
	 *
	 * @return The progress as a double between 0 and 1.
	 */
	double progress() const;
};



#endif /* INCLUDE_CONVOLVER_HPP_ */
