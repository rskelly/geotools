/*
 * convolve.hpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#ifndef INCLUDE_CONVOLVE_HPP_
#define INCLUDE_CONVOLVE_HPP_

#include <string>
#include <fstream>
#include <memory>

#include "util.hpp"
#include "reader.hpp"
#include "writer.hpp"

using namespace geo::util;
using namespace hlrg::reader;
using namespace hlrg::writer;

namespace hlrg {
namespace convolve {

/**
 * The kernel class holds the step values, plus the min
 * and max wavelengths for the range.
 */
class Kernel {
private:
	double m_wl;
	double m_fwhm;
	int m_window;
	int m_index;		///<! Index into the array of source wavelengths.
	std::vector<double> m_kernel;

public:

	/**
	 * Create a kernel around the given wavelength with the given
	 * full width at half maximum. The threshold gives the value
	 * of the Gaussian beyond which the kernel is not applied.
	 *
	 * \param size The number of kernel elements.
	 * \param bandDist The distance between the input bands.
	 * \param wl The centre wavelength.
	 * \param fwhm The full width at half maximum.
	 * \param window The size of the window; an odd integer.
	 */
	Kernel(double wl, double fwhm, int window, int index);

	int index() const;

	/**
	 * Set the wavelength.
	 *
	 * \param wl The wavelength.
	 */
	void setWavelength(double wl);

	/**
	 * Get the wavelength.
	 *
	 * \return The wavelength.
	 */
	double wl() const;

	/**
	 * Set the full width at half maximum.
	 *
	 * \param fwhm The full width at half maximum.
	 */
	void setFWHM(double fwhm);

	/**
	 * Get the full width at half maximum.
	 *
	 * \return The full width at half maximum.
	 */
	double fwhm() const;

	void setWindow(int window);

	int window() const;

	double apply(const std::vector<double>& intensities, const std::vector<double>& wavelengths, int idx) const;

	bool operator<(const Kernel& other) const;

	~Kernel();
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
	 * \param band The 1-based band index.
	 * \param wl The wavelength.
	 * \param fwhm The full width of half maximum.
	 */
	BandProp(int band, double wl, double fwhm);
};


/**
 * Represents a single band in the input/output spectra.
 */
class Band {
private:
	double m_wl;			///!< The wavelength.
	double m_scale;		 	///!< The scaled value of the intensity. This is used for calculations. Is identical to value by default.
	double m_shift;			///!< The amount to shift the band's wavelength designation by.
	int m_count;			///!< Tracks the number of accumulations; divide the value by this number.

public:

	/**
	 * Create a band.
	 *
	 * \param wl The wavelength.
	 * \param value The value or intensity.
	 */
	Band(double wl, double value);

	Band();

	/**
	 * Reset the counter.
	 */
	void reset();

	/**
	 * Return the count. The count is the number of times the value
	 * of the band has been updated by setValue.
	 *
	 * \return The count.
	 */
	int count() const;

	/**
	 * Set the scale factor.
	 *
	 * \param scale The scale factor.
	 */
	void setScale(double scale);

	/**
	 * Return the scale factor.
	 *
	 * \return The scale factor.
	 */
	double scale() const;

	/**
	 * Set the shift amount.
	 *
	 * \param shift The shift amount.
	 */
	void setShift(double shift);

	/**
	 * Get the shift amount.
	 *
	 * \return The shift amount.
	 */
	double shift() const;

	/**
	 * Return the Band's central wavelength.
	 *
	 * \return The Band's wavelength.
	 */
	double wl() const;

	/**
	 * Set he Band's wavelength.
	 *
	 * \param wl The wavelength.
	 */
	void setWl(double wl);

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
	char m_delim;									///<! The column delimiter.
	size_t m_count;									///<! The number of records in the file.
	int m_firstRow;									///<! The first row of data.
	int m_firstCol;									///<! The first column of data.
	int m_dateCol;									///<! The date column (or -1 if not used.)
	int m_timeCol;									///<! The timestamp column (or -1 if not used.)

	std::unique_ptr<GDALReader> m_raster;
	int m_col;
	int m_row;
	size_t m_rasterIdx;

	std::string m_projection;
	double m_trans[6];

	bool loadCSV(const std::string& filename, const std::string& delimiter);

	bool loadRaster(const std::string& filename, size_t memLimit);

public:
	std::vector<Band> bands;						///<! A list of the bands. This changes as the file is read through.
	std::vector<double> wavelengths;
	std::vector<double> intensities;
	std::map<std::string, std::string> properties;	///<! Properties read from the header block.
	std::string date;								///<! The date of the current row.
	long time;										///<! The timestamp of the current row.

	Spectrum();

	Spectrum(int firstRow, int firstCol, int dateCol, int timeCol);

	const std::string& projection() const;

	void transform(double* trans) const;

	/**
	 * Load the data file and read the header information.
	 * After load is called (and returns true), the next method must be
	 * called to populate the bands list with data.
	 *
	 * \param filename A data file.
	 * \param delimiter The column delimiter.
	 * \param firstRow The zero-based index of the first row of data.
	 * \param firstcol  The zero-based index of the first column of data.
	 * \return True if the file is loaded and has information in it.
	 */
	bool load(const std::string& filename, const std::string& delimiter, size_t memLimit);

	/**
	 * Return a guess of the size of the input file.
	 *
	 * \param filename A data file.
	 * \param delimiter The column delimiter.
	 */
	size_t inputSize(const std::string& filename, const std::string& delimiter);

	/**
	 * Return the number of spectra.
	 */
	size_t count() const;

	/**
	 * Advance the reader to the next line of data. This becomes the Spectrum's current state:
	 * the list of Bands contains data from the current row.
	 *
	 * \return True if a row has been read, false if there were none left.
	 */
	bool next();

	/**
	 * Sets up the given Spectrum with the same structure and
	 * wavelengths as the current one.
	 *
	 * \param spec A Spectrume.
	 */
	void setup(Spectrum& spec);

	/**
	 * Apply this scaling factor to every band's intensity value.
	 *
	 * \param scale A scale factor.
	 */
	void scale(double scale);

	/**
	 * Shift the band designations by this amount. This will be done
	 * before any other transformations.
	 *
	 * \param shift The amount to shift by, in wavelength units.
	 */
	void shift(double shift);

	/**
	 * Set the values of all bands to zero.
	 */
	void reset();

	/**
	 * Write the header of an output spectrum, between the given
	 * wavelengths.
	 *
	 * \param out An output stream.
	 * \param minWl The minimum wavelength to print.
	 * \param maxWl The maximum wavelength to print.
	 * \param delim The column delimiter.
	 */
	void writeHeader(std::ostream& out, double minWl, double maxWl, char delim);

	/**
	 * Write the bands of an output spectrum, between the given
	 * wavelengths.
	 *
	 * \param out An output stream.
	 * \param minWl The minimum wavelength to print.
	 * \param maxWl The maximum wavelength to print.
	 * \param delim The column delimiter.
	 */
	void write(std::ostream& out, double minWl, double maxWl, char delim);

	void write(GDALWriter* wtr, double minWl, double maxWl, int col, int row);

	int col() const;

	int row() const;

	std::unique_ptr<GDALReader>& raster();

};


/**
 * Read the convolution configuration file. The columns are
 * band #, wavelength, full width at half maximum.
 * Names are ignored; the order is important.
 */
class BandPropsReader {
private:
	std::map<int, BandProp> m_bandProps;	///<! The map of BandProp object.

public:
	double minWl;						///<! The minimum wavelength in the band definition file.
	double maxWl;						///<! The maximum wavelength in the band definition file.

	/**
	 * Load the band definition file.
	 *
	 * \param filename The filename.
	 * \param delimiter The column delimiter.
	 */
	void load(const std::string& filename, const std::string& delimiter);

	/**
	 * Returns a reference to the map containing the band configurations.
	 *
	 * \return A reference to the map containing the band configurations.
	 */
	const std::map<int, BandProp>& bands() const;

	/**
	 * Configure the given Spectrum with the bands represented
	 * by the band definition file.
	 *
	 * \param spec A Spectrum object.
	 */
	void configureSpectrum(Spectrum& spec);

};



/**
 * Forward declaration.
 */
class Convolve;

/**
 * An interface that provides implementors with the ability to receive
 * notifications from the Convolve.
 */
class ConvolveListener {
public:

	/**
	 * Called when the convolution has started.
	 *
	 * \param A Convolve.
	 */
	virtual void started(Convolve*) = 0;

	/**
	 * Called when the convolution status has updated.
	 *
	 * \param A Convolve.
	 */
	virtual void update(Convolve*) = 0;

	/**
	 * Called when the convolution has stopped.
	 *
	 * \param A Convolve.
	 */
	virtual void stopped(Convolve*) = 0;

	/**
	 * Called when the convolution has finished.
	 *
	 * \param A Convolve.
	 */
	virtual void finished(Convolve*) = 0;

	virtual ~ConvolveListener() {}
};


/**
 * Performs the convolution of a spectral file using Gaussian kernels
 * parameterized from a band definition file.
 */
class Convolve {
private:
	double m_progress;	///<! The current progress, between 0 and 1.

public:

	/**
	 * Run the convolve on the given files. The listener will receive callbacks.
	 *
	 * \param listener A ConvolveListener to receive updates.
	 * \param bandDef The band definition file.
	 * \param bandDefDelim The delimiter for the data file.
	 * \param spectra The spectral data file.
	 * \param spectraDelim The delimiter for the data file.
	 * \param spectraFirstRow The zero-based index of the first data row.
	 * \param spectraFirstCol The zero-based index of the first data column.
	 * \param spectraDateCol The zero-based index of a date string. -1 for none.
	 * \param spectraTimeCol The zero-based index of a timestamp. -1 for none.
	 * \param output The output file.
	 * \param outputDelim The delimiter for the data file.
	 * \param inputScale Scale every input spectral value by this much.
	 * \param tolerance A value that dictates how wide the Gaussian will be by providing a minimum threshold for the y-value.
	 * \param bandShift If given and non-zero, will cause the input band designations to be shifted by the given amount.
	 * \param memLimit The amount of memory that will trigger the use of file-backed storaged.
	 * \param threads The number of threads -- memLimit must be multiplied by this.
	 * \param running A reference to a boolean that is true so long as the processor should keep running.
	 */
	void run(ConvolveListener& listener,
			const std::string& bandDef, const std::string& bandDefDelim,
			const std::vector<std::string>& spectra, const std::string& spectraDelim,
			int spectraFirstRow, int spectraFirstCol, int spectraDateCol, int spectraTimeCol,
			const std::string& output, const std::string& outputDelim, FileType outputType,
			double inputScale, double tolerance, double bandShift, size_t memLimit, int threads, bool& running);

	/**
	 * Return a guess of the amount of memory needed to process the input file.
	 *
	 * \param spectra The spectral data file.
	 * \param spectraDelim The delimiter for the data file.
	 * \return The expected size of the processing dataset.
	 */
	size_t guess(const std::string& spectra, const std::string& spectraDelim);

	/**
	 * Return the progress as a double between 0 and 1.
	 *
	 * \return The progress as a double between 0 and 1.
	 */
	double progress() const;
};

} // convolve
} // hlrg

#endif /* INCLUDE_CONVOLVE_HPP_ */
