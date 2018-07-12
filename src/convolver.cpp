/*
 * convolver.cpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include "convolver.hpp"

#define PI 3.1415926535

/**
 * Compute the inverse of the Gaussian: retrieve the x that would give
 * the given y. X is returned as an absolute distance from x0.
 *
 * @param sigma The standard deviation of the function.
 * @param y The value of y for which to find x.
 * @param x0 The mean x.
 */
double invGaussian(double sigma, double y, double x0) {
	return sigma * std::sqrt(-2.0 * std::log(y));
}

/**
 * Compute the value of the Gaussian function for a given mean (x0) and
 * x, given sigma and magnitude 1.
 *
 * @param sigma The standard deviation of the function.
 * @param x The current x.
 * @param x0 The mean x.
 */
double gaussian(double sigma, double x, double x0) {
	return std::exp(-0.5 * std::pow((x - x0) / sigma, 2.0));
}

/**
 * The kernel class holds the step values, plus the min
 * and max wavelengths for the range.
 */
class Kernel {
public:
	double wl;
	double min;
	double max;
	std::vector<double> kernel;

	Kernel(double wl, double min, double max) :
		wl(wl), min(min), max(max) {}

	Kernel() : Kernel(0, 0, 0) {}

	void resize(int size) {
		kernel.resize(size);
		kernel[0] = 0;
		kernel[size - 1] = 0;
	}

	double& operator[](size_t idx) {
		return kernel[idx];
	}

	/**
	 * Normalize the kernal so that its elements sum to 1.
	 */
	void normalize() {
		double s = 0;
		for(double d : kernel)
			s += d;
		for(double& d : kernel)
			d /= s;
	}

	/**
	 * Build a 1-dimensional kernel contaning the discretized Gaussian with
	 * size elements. The Gaussian is parameterized by its mean (the bandwidth) and
	 * standard deviation, which is computed from the full width at half magnitude (FWHM)
	 * using the function 2*sqrt(2*ln(2)) * sigma (https://en.wikipedia.org/wiki/Full_width_at_half_maximum).
	 *
	 * The kernel is normalized to sum to 1.
	 *
	 * @param size The number of elements in the kernel (one element added if even, so that there's a zero at each end and the peak in the middle).
	 * @param wl The wavelength in whatever units; must be uniform for all sources.
	 * @param fwhm Full width at half maximum.
	 * @param threshold The minimum value of the function, which will delineate the "ends" at index 0, and size.
	 */
	void build(int size, double wl, double fwhm, double threshold) {
		if(size % 2 == 0) ++size;
		// Half-width at half maximum.
		double hwhm = fwhm * 0.5;
		// Recover the std. dev from the hwhm.
		double sigma = hwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
		// Get the distance from wl for the threshold value.
		double inv = invGaussian(sigma, threshold, wl);
		// Find the bounds for the curve, and the step width.
		this->wl = wl;
		min = wl - inv;
		max = wl + inv;
		double step = (max - min) / (size - 1);
		// Compute the steps.
		kernel.resize(size);
		for(int i = 0; i < size; ++i)
			kernel[i] = gaussian(sigma, min + step * i, wl);
		// Normalize the kernel.
		normalize();
	}
};


/**
 * Contains the information about each band that will be used to build
 * the convolution kernel.
 */
class BandProp {
public:
	int band;
	double wl;
	double fwhm;
	BandProp(int band, double wl, double fwhm) :
		band(band), wl(wl), fwhm(fwhm) {}
};

class Band {
public:
	double wl;
	double value;
	Band(double wl, double value) :
		wl(wl), value(value) {}
	Band() : Band(0, 0) {}
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
	std::ifstream m_input;
	std::string m_buf;
	size_t m_count;

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
	bool load(const std::string& filename) {
		bands.clear();
		m_input.open(filename, std::ios::in);
		std::string part;
		// Try to skip the header. If it doesn't work quit.
		if(!std::getline(m_input, m_buf))
			return false;
		// Run over the rows.
		bool header = false; // The band wl header hasn't been read yet.
		while(std::getline(m_input, m_buf)) {
			//count = (pos = buf.find(',')) == std::string::npos ? buf.size() : buf.size() - pos;
			//fwhm = atof(buf.substr(0, count).c_str()); // May or may not be a comma at the end.
			if(!header && m_buf.find(':') < std::string::npos) {
				// Process the colon-delimited headers for properties.
				std::stringstream ss(m_buf);
				std::getline(ss, part, ':');
				std::getline(ss, m_buf);
				properties[part] = m_buf;
			} else if(m_buf.find(">>>") < std::string::npos) {
				// Skip the divider.
				continue;
			} else if(!header) {
				// Parse the wavelengths out of the header section. The first two columns (date, time) are empty.
				std::stringstream ss(m_buf);
				std::getline(ss, part, '\t');
				std::getline(ss, part, '\t');
				while(std::getline(ss, part, '\t'))
					bands.emplace_back(atof(part.c_str()), 0);
				header = true;
			} else {
				// We're at the start of data.
				size_t pos = m_input.tellg();
				std::string dummy;
				m_count = 1;
				while(std::getline(m_input, dummy))
					++m_count;
				m_input.close();
				m_input.open(filename, std::ios::in);
				if(!m_input.seekg(pos).good())
					return false;
				break;
			}
		}
		return true;
	}

	/**
	 * Return the number of spectra.
	 */
	size_t count() const {
		return m_count;
	}

	/**
	 * Advance the reader to the next line of data. This becomes the Spectrum's current state:
	 * the list of Bands contains data from the current row.
	 *
	 * @return True if a row has been read, false if there were none left.
	 */
	bool next() {
		if(m_buf.empty())
			return false;

		std::stringstream ss(m_buf);
		std::string part;
		std::getline(ss, date, '\t');
		std::getline(ss, part, '\t');
		time = std::strtol(part.c_str(), nullptr, 10);
		size_t i = 0;
		while(std::getline(ss, part, '\t'))
			bands[i++].value = std::strtod(part.c_str(), nullptr);

		std::getline(m_input, m_buf);
		return true;
	}

	void setup(Spectrum& spec) {
		// Instantiate the bands on the new spectrum.
		for(const Band& b : bands)
			spec.bands.emplace_back(b.wl, 0);
	}

	void convolve(Kernel& kernel, Band& band) {
		// Find the index nearest the centre of the kernel.
		size_t idx = 0;
		for(size_t i = 0; i < bands.size(); ++i) {
			if(kernel.min < bands[i].wl) break;
			idx = i;
		}
		// Convolve the bands.
		for(size_t i = 0; i < kernel.kernel.size(); ++i)
			band.value += bands[i + idx].value * kernel[i];
	}

	void reset() {
		for(Band& b : bands)
			b.value = 0;
	}

	void writeHeader(std::ostream& out, double minWl, double maxWl) {
		out << "date,timestamp";
		for(const Band& b : bands) {
			if(b.wl >= minWl && b.wl <= maxWl)
				out << "," << b.wl;
		}
		out << "\n";
	}

	void write(std::ostream& out, double minWl, double maxWl) {
		out << date << "," << time;
		for(const Band& b : bands) {
			if(b.wl >= minWl && b.wl <= maxWl)
				out << "," << b.value;
		}
		out << "\n";
	}

};

/**
 * Read the convolution configuration file. The columns are
 * band #, wavelength, full width at half maximum.
 * Names are ignored; the order is important.
 */
class BandPropsReader {
public:
	std::map<int, BandProp> bandProps;
	double minWl;
	double maxWl;

	void load(const std::string& filename) {
		minWl = std::numeric_limits<double>::max();
		maxWl = std::numeric_limits<double>::lowest();

		std::ifstream input(filename, std::ios::in);
		std::string buf, part;
		// Try to skip the header. If it doesn't work quit.
		if(!std::getline(input, buf))
			return;
		// Run over the rows.
		int band;
		double wl, fwhm;
		size_t pos, count;
		while(std::getline(input, buf)) {
			std::stringstream ss(buf);
			std::getline(ss, buf, ',');
			band = std::strtol(buf.c_str(), nullptr, 10);
			std::getline(ss, buf, ',');
			wl = std::strtod(buf.c_str(), nullptr);
			std::getline(ss, buf);
			count = (pos = buf.find(',')) == std::string::npos ? buf.size() : buf.size() - pos;
			fwhm = std::strtod(buf.substr(0, count).c_str(), nullptr); // May or may not be a comma at the end.
			bandProps.emplace(std::piecewise_construct, std::forward_as_tuple(band), std::forward_as_tuple(band, wl, fwhm));
			if(wl < minWl) minWl = wl;
			if(wl > maxWl) maxWl = wl;
		}
	}

	std::vector<int> bands() const {
		std::vector<int> bands;
		for(const auto& p : bandProps)
			bands.push_back(p.first);
		return bands;
	}

	void configureKernel(Kernel& kernel, int band, int size, double threshold) {
		if(bandProps.find(band) == bandProps.end())
			throw std::runtime_error("Band not found: " + std::to_string(band));
		const BandProp& p = bandProps.at(band);
		kernel.build(size, p.wl, p.fwhm, threshold);
	}

	void configureSpectrum(Spectrum& spec) {
		spec.bands.resize(bandProps.size());
		size_t i = 0;
		for(const auto& p : bandProps)
			spec.bands[i++].wl = p.second.wl;
	}
};

void Convolver::run(ConvolverListener* listener, const std::string& bandDef, const std::string& spectra, const std::string& output) {
	listener->started(this);

	Kernel kernel;
	BandPropsReader rdr;
	rdr.load(bandDef);

	Spectrum spec;
	spec.load(spectra);

	Spectrum out;
	rdr.configureSpectrum(out);

	std::ofstream outstr(output, std::ios::out);
	size_t complete = 0, count = spec.count();
	bool header = false;
	while(spec.next()) {
		out.reset();
		out.date = spec.date;
		out.time = spec.time;
		for(int b : rdr.bands()) {
			rdr.configureKernel(kernel, b, 40, 0.001);
			spec.convolve(kernel, out.bands[b - 1]);
		}
		if(!header) {
			out.writeHeader(outstr, rdr.minWl, rdr.maxWl);
			header = true;
		}
		out.write(outstr, rdr.minWl, rdr.maxWl);
		m_progress = (double) complete++ / count;
		listener->update(this);
	}
	listener->stopped(this);
}

void Convolver::cancel() {

}

double Convolver::progress() const {
	return m_progress;
}

