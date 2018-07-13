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
#include <limits>

#include "convolver.hpp"

constexpr double PI = 3.1415926535;

/**
 * Compute the inverse of the Gaussian: retrieve the x that would give
 * the given y. X is returned as an absolute distance from x0.
 *
 * @param sigma The standard deviation of the function.
 * @param y The value of y for which to find x.
 * @param x0 The mean x.
 */
double invGaussian(double sigma, double y, double x0) {
	return sigma * std::sqrt(-2.0 * std::log(y * 2.0 * PI)) + x0;
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
	return  1.0 / std::sqrt(2.0 * PI) * std::exp(-0.5 * std::pow((x - x0) / sigma, 2.0));
}


Kernel::Kernel(double wl, double fwhm, double threshold) :
	m_wl(wl),
	m_fwhm(fwhm),
	m_threshold(threshold) {
	init();
}

void Kernel::setWavelength(double wl) {
	m_wl = wl;
	init();
}

double Kernel::wl() const {
	return m_wl;
}

void Kernel::setFWHM(double fwhm) {
	m_fwhm = fwhm;
	init();
}

double Kernel::fwhm() const {
	return m_fwhm;
}

void Kernel::init() {
	m_sigma = ((0.5 * m_fwhm) / (2.0 * std::sqrt(2.0 * std::log(2.0))));
	double w = invGaussian(m_sigma, m_threshold, m_wl);
	m_halfWidth = std::abs(m_wl - w);
}

double Kernel::halfWidth() const {
	return m_halfWidth;
}

double Kernel::operator()(double wl0) const {
	return gaussian(m_sigma, wl0, m_wl);
}

double Kernel::threshold() const {
	return m_threshold;
}



BandProp::BandProp(int band, double wl, double fwhm) :
		band(band), wl(wl), fwhm(fwhm) {
}



Band::Band(double wl, double value) :
			wl(wl), value(value), scaledValue(value) {
}

Band::Band() : Band(0, 0) {
}



bool Spectrum::load(const std::string& filename) {
	// Clear any existing bands list.
	bands.clear();
	// Open the input file for reading.
	m_input.open(filename, std::ios::in);

	// Try to skip the header. If it doesn't work quit.
	if(!std::getline(m_input, m_buf))
		return false;

	std::string buf;		// A temporary buffer.
	bool header = false; 	// False when the band wl header hasn't been read yet.

	// Run over the rows.
	while(std::getline(m_input, m_buf)) {
		if(!header && m_buf.find(':') < std::string::npos) {
			// Process the colon-delimited headers for properties.
			std::stringstream ss(m_buf);
			std::getline(ss, buf, ':');
			std::getline(ss, m_buf);
			properties[buf] = m_buf;
		} else if(m_buf.find(">>>") < std::string::npos) {
			// Skip the divider.
			continue;
		} else if(!header) {
			// Parse the wavelengths out of the header section. The first two columns (date, time) are empty.
			std::stringstream ss(m_buf);
			std::getline(ss, buf, '\t');
			std::getline(ss, buf, '\t');
			while(std::getline(ss, buf, '\t'))
				bands.emplace_back(atof(buf.c_str()), 0);
			header = true;
		} else {
			// We're at the start of data. Get the position, and count the remaining rows.
			size_t pos = m_input.tellg();
			std::string dummy;
			m_count = 1;
			while(std::getline(m_input, dummy))
				++m_count;

			// Close and open the file, and seek back to where we left off.
			// The first row is already in the buffer.
			m_input.close();
			m_input.open(filename, std::ios::in);
			if(!m_input.seekg(pos).good())
				return false;
			break;
		}
	}
	return true;
}

size_t Spectrum::count() const {
	return m_count;
}

bool Spectrum::next() {
	// No buffer was read on the previous read. We're done.
	if(m_buf.empty())
		return false;
	// Read the date, timestamp and values out of the line.
	std::stringstream ss(m_buf);
	std::string part;
	std::getline(ss, date, '\t');
	std::getline(ss, part, '\t');
	time = std::strtol(part.c_str(), nullptr, 10);
	size_t i = 0;
	while(std::getline(ss, part, '\t'))
		bands[i++].value = std::strtod(part.c_str(), nullptr);
	// Read the next buffer. If it fails, we'll find out on the next call to next.
	std::getline(m_input, m_buf);
	return true;
}

void Spectrum::setup(Spectrum& spec) {
	// Instantiate the bands on the new spectrum.
	for(const Band& b : bands)
		spec.bands.emplace_back(b.wl, 0);
}

void Spectrum::convolve(Kernel& kernel, Band& band) {
	double min = band.wl - kernel.halfWidth();
	double max = band.wl + kernel.halfWidth();
	// Find the index nearest the centre of the kernel.
	size_t idx0 = 0;
	size_t idx1 = bands.size();
	for(size_t i = 0; i < bands.size(); ++i) {
		if(min > bands[i].wl)
			idx0 = i;
		else break;
	}
	for(size_t i = idx0; i < bands.size(); ++i) {
		if(max < bands[i].wl) {
			idx1 = i;
			break;
		}
	}

	// Convolve the bands.
	for(size_t i = idx0; i < idx1; ++i) {
		double v = kernel(bands[i].wl);
		band.value += bands[i].scaledValue * v;
	}
}

void Spectrum::reset() {
	for(Band& b : bands)
		b.scaledValue = b.value = 0;

}

void Spectrum::writeHeader(std::ostream& out, double minWl, double maxWl) {
	out << "date,timestamp";
	for(const Band& b : bands) {
		if(b.wl >= minWl && b.wl <= maxWl)
			out << "," << b.wl;
	}
	out << "\n";
}

void Spectrum::write(std::ostream& out, double minWl, double maxWl) {
	out << date << "," << time;
	for(const Band& b : bands) {
		if(b.wl >= minWl && b.wl <= maxWl)
			out << "," << b.value;
	}
	out << "\n";
}

void Spectrum::scale(double scale) {
	for(Band& b : bands)
		b.scaledValue = b.value * scale;
}

void BandPropsReader::load(const std::string& filename) {
	// Set the initial limits. These will be refined.
	minWl = std::numeric_limits<double>::max();
	maxWl = std::numeric_limits<double>::lowest();
	// Open file, initialize buffers.
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

std::vector<int> BandPropsReader::bands() const {
	std::vector<int> bands;
	for(const auto& p : bandProps)
		bands.push_back(p.first);
	return bands;
}

void BandPropsReader::configureKernel(Kernel& kernel, int band) {
	// If the band is not found, raise an error.
	if(bandProps.find(band) == bandProps.end())
		throw std::runtime_error("Band not found: " + std::to_string(band));
	// Get the band and set the properties.
	const BandProp& p = bandProps.at(band);
	kernel.setWavelength(p.wl);
	kernel.setFWHM(p.fwhm);
}

void BandPropsReader::configureSpectrum(Spectrum& spec) {
	// Resize the spectrum to the number of bands in the configuration.
	spec.bands.resize(bandProps.size());
	// Set all the wavelengths on the bands.
	size_t i = 0;
	for(const auto& p : bandProps)
		spec.bands[i++].wl = p.second.wl;
}



void Convolver::run(ConvolverListener& listener,
		const std::string& bandDef, const std::string& spectra, const std::string& output,
		double inputScale, bool& running) {

	// Notify a listener.
	listener.started(this);

	// Configure the kernel; only the threshold value available at this time.
	Kernel kernel(0, 0, 0.0001);

	// Load the band properties.
	BandPropsReader rdr;
	rdr.load(bandDef);

	// Load the spectrum.
	Spectrum spec;
	spec.load(spectra);
	spec.scale(inputScale);

	// Configure the output (convolved) spectrum.
	Spectrum out;
	rdr.configureSpectrum(out);

	// Open the output file.
	std::ofstream outstr(output, std::ios::out);

	// Run the convolution record-by-record.
	size_t complete = 0, count = spec.count(); // Status counters.
	bool header = false;
	while(running && spec.next()) {
		out.reset();
		out.date = spec.date;
		out.time = spec.time;
		for(int b : rdr.bands()) {
			// For each band definition, configure the kernel.
			rdr.configureKernel(kernel, b);
			// Convolve the spectrum onto the given band.
			spec.convolve(kernel, out.bands[b - 1]);
			if(!running) break;
		}
		if(!header) {
			// Write the header if this is the first iteration.
			out.writeHeader(outstr, rdr.minWl, rdr.maxWl);
			header = true;
		}
		// Write the record.
		out.write(outstr, rdr.minWl, rdr.maxWl);
		// Update progress.
		m_progress = (double) complete++ / count;
		listener.update(this);
		if(!running) break;
	}
	// Notify listener of completion.
	listener.finished(this);
}

double Convolver::progress() const {
	return m_progress;
}

