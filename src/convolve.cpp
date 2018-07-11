/*
 * convolve.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: rob
 */

#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

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
	double min;
	double max;
	std::vector<double> kernel;

	Kernel(double min, double max) :
		min(min), max(max) {}

	Kernel() : Kernel(0, 0) {}

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
			band = atoi(buf.c_str());
			std::getline(ss, buf, ',');
			wl = atof(buf.c_str());
			std::getline(ss, buf);
			count = (pos = buf.find(',')) == std::string::npos ? buf.size() : buf.size() - pos;
			fwhm = atof(buf.substr(0, count).c_str()); // May or may not be a comma at the end.
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

};


class Band {
public:
	double wl;
	double value;
	Band(double wl, double value) :
		wl(wl), value(value) {}
	Band() : Band(0, 0) {}
};

class Spectrum {
private:
	std::ifstream m_input;
	std::string m_buf;

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
				break;
			}
		}
		return true;
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
		time = atoi(part.c_str());
		size_t i = 0;
		while(std::getline(ss, part, '\t'))
			bands[i++].value = atof(part.c_str());

		std::getline(m_input, m_buf);
		return true;
	}

	void setup(Spectrum& spec) {
		// Instantiate the bands on the new spectrum.
		for(const Band& b : bands)
			spec.bands.emplace_back(b.wl, 0);
	}

	void convolve(Kernel& kernel, Spectrum& spec) {
		// Convolve the bands.
		for(size_t i = 0; i < bands.size(); ++i) {
			const Band& b = bands[i];
			if(b.wl >= kernel.min && b.wl <= kernel.max) {
				Band& b0 = spec.bands[i];
				for(double k : kernel.kernel)
					b0.value += b.value * k;
			}
		}
	}

	void reset() {
		for(Band& b : bands)
			b.value = 0;
	}

	void writeHeader(std::ostream& out, double minWl, double maxWl) {
		for(const Band& b : bands) {
			if(b.wl >= minWl && b.wl <= maxWl)
				out << b.wl << ",";
		}
		out << "\n";
	}

	void write(std::ostream& out, double minWl, double maxWl) {
		for(const Band& b : bands) {
			if(b.wl >= minWl && b.wl <= maxWl)
				out << b.value << ",";
		}
		out << "\n";
	}

};

int main(int argc, char** argv) {

	Kernel kernel;
	BandPropsReader rdr;
	rdr.load("../data/Headwall_NANO_Bands.csv");

	Spectrum spec;
	Spectrum output;
	spec.load("../data/AbsoluteIrradiance_10-35-48-938.txt");
	spec.setup(output);

	bool header = false;
	int i = 0;
	while(spec.next()) {
		output.reset();
		for(int b : rdr.bands()) {
			rdr.configureKernel(kernel, b, 40, 0.001);
			spec.convolve(kernel, output);
		}
		if(!header) {
			output.writeHeader(std::cout, rdr.minWl, rdr.maxWl);
			header = true;
		}
		output.write(std::cout, rdr.minWl, rdr.maxWl);
		if(++i > 100)
			break;
	}

	/*
	int size = 40;
	Kernel kernel;
	kernel.build(size, 100, 5, 0.001);

	Spectrum orig;
	orig.bands.emplace_back(98, 0.1);
	orig.bands.emplace_back(99, 0.2);
	orig.bands.emplace_back(100, 0.5);
	orig.bands.emplace_back(101, 0.3);

	Spectrum conv;
	orig.convolve(kernel, conv);

	std::cerr << kernel.min << ", " << kernel.max << "\n";
	for(Band& b : conv.bands)
		std::cerr << b.wl << ", " << b.value << "\n";
	*/

	return 0;
}
