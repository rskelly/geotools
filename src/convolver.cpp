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

using namespace hlrg;

namespace {

	/**
	 * Compute the inverse of the Gaussian: retrieve the x that would give
	 * the given y. X is returned as an absolute distance from x0.
	 *
	 * @param sigma The standard deviation of the function.
	 * @param y The value of y for which to find x.
	 * @param x0 The mean x.
	 */
	inline double invGaussian(double sigma, double y, double x0) {
		return sigma * std::sqrt(-2.0 * std::log(y * sigma * std::sqrt(2.0 * M_PI))) + x0;
	}

	/**
	 * Compute the value of the Gaussian function for a given mean (x0) and
	 * x, given sigma and magnitude 1.
	 *
	 * @param sigma The standard deviation of the function.
	 * @param x The current x.
	 * @param x0 The mean (expected) x.
	 */
	inline double gaussian(double sigma, double x, double x0) {
		return  1.0 / (sigma * std::sqrt(2.0 * M_PI)) * std::exp(-0.5 * std::pow((x - x0) / sigma, 2.0));
	}

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
	// Derive the std dev from the FWHM.
	m_sigma = m_fwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
	// Get the distance of the threshold on x from threshold y.
	double w = invGaussian(m_sigma, m_threshold, m_wl);
	// The half-width of the curve at the threshold.
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
		m_wl(wl),
		m_value(value),
		m_scale(1),
		m_shift(0),
		m_count(0) {
}

Band::Band() : Band(0, 0) {
}

void Band::setValue(double value) {
	m_value = value;
	++m_count;
}

double Band::value() const {
	return m_value;
}

void Band::setShift(double shift) {
	m_shift = shift;
}

double Band::shift() const {
	return m_shift;
}

double Band::normalizedValue() const {
	return m_count > 0 ? m_value / m_count : std::numeric_limits<double>::quiet_NaN();
}

double Band::scaledValue() const {
	return m_value * m_scale;
}

void Band::setScale(double scale) {
	m_scale = scale;
}

double Band::scale() const {
	return m_scale;
}

double Band::wl() const {
	return m_wl + m_shift;
}

void Band::setWl(double wl) {
	m_wl = wl;
}

void Band::reset() {
	m_count = 0;
}

int Band::count() const {
	return m_count;
}

void _stripws(std::string& buf) {
	size_t j = 0;
	for(size_t i = 0; i < buf.size(); ++i) {
		char c = buf[i];
		if(c != '\n' && c != '\r' && c != ' ')
			buf[j++] = c;
	}
	buf.resize(j);
}

void _stripcr(std::string& buf) {
	size_t j = 0;
	for(size_t i = 0; i < buf.size(); ++i) {
		char c = buf[i];
		if(c != '\r') buf[j++] = c;
	}
	buf.resize(j);
}

Spectrum::Spectrum(int firstRow, int firstCol, int dateCol, int timeCol) :
		m_delim(','),
		m_count(0),
		m_firstRow(firstRow), m_firstCol(firstCol),
		m_dateCol(dateCol), m_timeCol(timeCol),
		time(0) {}

Spectrum::Spectrum() : Spectrum(0, 0, -1, -1) {}

bool Spectrum::load(const std::string& filename, const std::string& delimiter) {
	// Clear any existing bands list.
	bands.clear();
	// Open the input file for reading.
	m_input.open(filename, std::ios::in);
	// Set the delimiter
	m_delim = delimiter[0];

	// Try to skip extraneous rows.
	for(int r = 0; r < m_firstRow; ++r) {
		if(!std::getline(m_input, m_buf))
			return false;
	}

	std::string buf;		// A temporary buffer.
	bool header = false; 	// False when the band wl header hasn't been read yet.

	// Run over the rows.
	while(std::getline(m_input, m_buf)) {
		_stripcr(m_buf); // Strip carriage return.
		if(m_buf.size() == 0) {
			// If the line contains no information, skip it.
			continue;
		} else if(!header) {
			if(std::string::npos == m_buf.find(m_delim))
				throw std::runtime_error("The selected delimiter wasn't found in this file.");
			// Parse the wavelengths out of the header section. The first two columns (date, time) are empty.
			// Skip the unneeded columns.
			std::stringstream ss(m_buf);
			for(int c = 0; c < m_firstCol; ++c)
				std::getline(ss, buf, m_delim);
			while(std::getline(ss, buf, m_delim))
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

	// Read the date string.
	{
		std::stringstream ss(m_buf);
		for(int c = 0; c <= m_dateCol; ++c)
			std::getline(ss, date, m_delim);
	}

	// Read the timestamp.
	{
		std::stringstream ss(m_buf);
		std::string part;
		for(int c = 0; c <= m_timeCol; ++c)
			std::getline(ss, part, m_delim);
		time = std::strtol(part.c_str(), nullptr, 10);
	}

	// Read the data.
	{
		size_t i = 0;
		std::stringstream ss(m_buf);
		std::string part;
		for(int c = 0; c < m_firstCol; ++c)
			std::getline(ss, part, m_delim);
		while(std::getline(ss, part, m_delim))
			bands[i++].setValue(std::strtod(part.c_str(), nullptr));
	}

	// Read the next buffer. If it fails, we'll find out on the next call to next.
	if(!m_input.eof() && m_input.good()) {
		std::getline(m_input, m_buf);
		if(!m_buf.empty() && std::string::npos == m_buf.find(m_delim))
			throw std::runtime_error("The selected delimiter wasn't found in this file.");
	} else {
		m_buf = "";
	}
	return true;
}

void Spectrum::setup(Spectrum& spec) {
	// Instantiate the bands on the new spectrum.
	for(const Band& b : bands)
		spec.bands.emplace_back(b.wl(), 0);
}

void Spectrum::convolve(Kernel& kernel, Band& band) {
	// First, figure out the start and end indices of the input, given
	// the width of the kernel function out to the threshold.
	double min = band.wl() - kernel.halfWidth() * 2;
	double max = band.wl() + kernel.halfWidth() * 2;
	size_t idx0 = 0;
	size_t idx1 = bands.size();
	for(size_t i = 0; i < bands.size(); ++i) {
		if(min > bands[i].wl()) {
			idx0 = i;
		} else {
			break;
		}
	}
	for(size_t i = idx0; i < bands.size(); ++i) {
		if(max < bands[i].wl()) {
			idx1 = i;
			break;
		}
	}
	// Build a temporary "kernel" to store the computed values from the
	// kernel function.
	std::vector<double> k(idx1 - idx0 + 1);
	double sum = 0;
	for(size_t i = idx0; i <= idx1; ++i)
		sum += k[i - idx0] = kernel(bands[i].wl());
	// Normalize the kernel values.
	for(size_t i = idx0; i <= idx1; ++i)
		k[i - idx0] /= sum;
	// Convolve the bands using the temporary kernel.
	for(size_t i = idx0; i <= idx1; ++i) {
		double v = k[i - idx0];
		double w = bands[i].scaledValue();
		band.setValue(band.value() + v * w);
	}
}

void Spectrum::reset() {
	for(Band& b : bands) {
		b.setValue(0);
		b.reset();
	}
}

void Spectrum::writeHeader(std::ostream& out, double minWl, double maxWl, char delim) {
	out << "date,timestamp";
	for(const Band& b : bands) {
		if(b.wl() >= minWl && b.wl() <= maxWl)
			out << delim << b.wl();
	}
	out << "\n";
}

void Spectrum::write(std::ostream& out, double minWl, double maxWl, char delim) {
	out << date << delim << time;
	for(const Band& b : bands) {
		if(b.wl() >= minWl && b.wl() <= maxWl)
			out << delim << b.value();
	}
	out << "\n";
}

void Spectrum::scale(double scale) {
	for(Band& b : bands)
		b.setScale(scale);
}

void Spectrum::shift(double shift) {
	for(Band& b : bands)
		b.setShift(shift);
}

void BandPropsReader::load(const std::string& filename, const std::string& delimiter) {
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
	char delim = delimiter[0];
	size_t pos, count;
	while(std::getline(input, buf)) {
		std::stringstream ss(buf);
		std::getline(ss, buf, delim);
		band = std::strtol(buf.c_str(), nullptr, 10);
		std::getline(ss, buf, delim);
		wl = std::strtod(buf.c_str(), nullptr);
		std::getline(ss, buf);
		count = (pos = buf.find(delim)) == std::string::npos ? buf.size() : pos;
		fwhm = std::strtod(buf.substr(0, count).c_str(), nullptr); // May or may not be a comma at the end.
		bandProps.emplace(std::piecewise_construct, std::forward_as_tuple(band), std::forward_as_tuple(band, wl, fwhm));
		if(wl < minWl) minWl = wl;
		if(wl > maxWl) maxWl = wl;
	}
	m_bands.clear();
	for(const auto& it : bandProps)
		m_bands.push_back(it.first);
}

const std::vector<int>& BandPropsReader::bands() const {
	return m_bands;
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
		spec.bands[i++].setWl(p.second.wl);
}

void Convolver::run(ConvolverListener& listener,
		const std::string& bandDef,
		const std::string& bandDefDelim,
		const std::string& spectra,
		const std::string& spectraDelim,
		int spectraFirstRow, int spectraFirstCol,
		int spectraDateCol, int spectraTimeCol,
		const std::string& output,
		const std::string& outputDelim,
		double inputScale, double tolerance, double bandShift, bool& running) {

	// Notify a listener.
	listener.started(this);

	// Configure the kernel.
	Kernel kernel(0, 0, tolerance);

	// Load the band properties.
	BandPropsReader rdr;
	rdr.load(bandDef, bandDefDelim);

	// Load the spectrum.
	Spectrum spec(spectraFirstRow, spectraFirstCol, spectraDateCol, spectraTimeCol);
	spec.load(spectra, spectraDelim);
	spec.shift(bandShift);
	spec.scale(inputScale);

	// Configure the output (convolved) spectrum.
	Spectrum out;
	rdr.configureSpectrum(out);

	// Open the output file.
	std::ofstream outstr(output, std::ios::out|std::ios::trunc);

	// Run the convolution record-by-record.
	size_t complete = 0, count = spec.count(); // Status counters.
	bool header = false;
	char delim = outputDelim[0];
	const std::vector<int>& bands = rdr.bands();
	while(running && spec.next()) {
		out.reset();
		out.date = spec.date;
		out.time = spec.time;
		for(size_t i = 0; i < bands.size(); ++i) {
			// For each band definition, configure the kernel.
			rdr.configureKernel(kernel, bands[i]);
			// Convolve the spectrum onto the given band.
			spec.convolve(kernel, out.bands[i]);
			if(!running) break;
		}
		if(!header) {
			// Write the header if this is the first iteration.
			out.writeHeader(outstr, rdr.minWl, rdr.maxWl, delim);
			header = true;
		}
		// Write the record.
		out.write(outstr, rdr.minWl, rdr.maxWl, delim);
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

