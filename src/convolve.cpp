/*
 * convolve.cpp
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

#include "convolve.hpp"
#include "writer.hpp"

using namespace hlrg::convolve;
using namespace hlrg::writer;

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


	/*
	void _stripws(std::string& buf) {
		size_t j = 0;
		for(size_t i = 0; i < buf.size(); ++i) {
			char c = buf[i];
			if(c != '\n' && c != '\r' && c != ' ')
				buf[j++] = c;
		}
		buf.resize(j);
	}
	 */

	void _stripcr(std::string& buf) {
		size_t j = 0;
		for(size_t i = 0; i < buf.size(); ++i) {
			char c = buf[i];
			if(c != '\r') buf[j++] = c;
		}
		buf.resize(j);
	}

}

Kernel::Kernel(double wl, double fwhm, int window, int index) :
	m_wl(wl),
	m_fwhm(fwhm),
	m_window(window),
	m_index(index) {

	m_sigma = std::sqrt(fwhm / (2.0 * std::sqrt(2.0 * std::log(2.0))));
	m_norm = 1.0 / (m_sigma * std::sqrt(2.0 * M_PI));
}

int Kernel::index() const {
	return m_index;
}

void Kernel::setWavelength(double wl) {
	m_wl = wl;
}

double Kernel::wl() const {
	return m_wl;
}

void Kernel::setFWHM(double fwhm) {
	m_fwhm = fwhm;
}

double Kernel::fwhm() const {
	return m_fwhm;
}

void Kernel::setWindow(int window) {
	m_window = window;
}

int Kernel::window() const {
	return m_window;
}

double Kernel::apply(const std::vector<double>& intensities, const std::vector<double>& wavelengths, int idx) const {
	double out = 0;
	for(int i = std::max(0, idx - m_window / 2); i < std::min((int) intensities.size(), idx + m_window / 2 + 1); ++i)
		out += m_norm * std::exp(-0.5 * std::pow((wavelengths[i] - m_wl) / m_sigma, 2.0));
	return out;
}

bool Kernel::operator<(const Kernel& other) const {
	return wl() < other.wl();
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


Spectrum::Spectrum(int firstRow, int firstCol, int dateCol, int timeCol) :
		m_delim(','),
		m_count(0),
		m_firstRow(firstRow), m_firstCol(firstCol),
		m_dateCol(dateCol), m_timeCol(timeCol),
		m_col(0), m_row(0),
		m_rasterIdx(0),
		time(0) {}

Spectrum::Spectrum() : Spectrum(0, 0, -1, -1) {}

size_t Spectrum::inputSize(const std::string& filename, const std::string& delimiter) {
	std::ifstream input(filename);
	input.seekg(0, std::ios::end);
	return input.tellg();
}

void Spectrum::update() {
	for(size_t i = 0; i < intensities.size(); ++i)
		bands[i].setValue(intensities[i]);
}

bool Spectrum::load(const std::string& filename, const std::string& delimiter, size_t memLimit) {
	// Clear any existing bands list.
	bands.clear();
	if(getFileType(filename) == FileType::CSV) {
		return loadCSV(filename, delimiter);
	} else {
		return loadRaster(filename, memLimit);
	}
}

bool Spectrum::loadRaster(const std::string& filename, size_t memLimit) {

	m_raster.reset(new GDALReader(filename, memLimit));
	m_projection = m_raster->projection();
	m_raster->transform(m_trans);
	m_raster->remap();
	std::vector<double> bandRange = m_raster->getWavelengths();

	for(double b : bandRange) {
		bands.emplace_back(b, 0);
		wavelengths.push_back(b);
	}

	intensities.resize(wavelengths.size());

	m_count = m_raster->cols() * m_raster->rows();

	return true;
}

bool Spectrum::loadCSV(const std::string& filename, const std::string& delimiter) {

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
			while(std::getline(ss, buf, m_delim)) {
				double b = std::stod(buf.c_str(), 0);
				bands.emplace_back(b, 0);
				wavelengths.push_back(b);
			}

			intensities.resize(wavelengths.size());

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

int Spectrum::col() const {
	return m_col;
}

int Spectrum::row() const {
	return m_row;
}

bool Spectrum::next() {

	if(m_raster.get()) {

		if(m_rasterIdx < m_count) {
			m_col = m_rasterIdx % m_raster->cols();
			m_row = m_rasterIdx / m_raster->cols();
			std::vector<double> values;
			m_raster->mapped(m_col, m_row, values);

			for(size_t i = 0; i < bands.size(); ++i) {
				bands[i].setValue(values[i]);
				intensities[i] = values[i];
			}

			++m_rasterIdx;
			return true;
		} else {
			return false;
		}

	} else {

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
			while(std::getline(ss, part, m_delim)) {
				double v = std::strtod(part.c_str(), nullptr);
				bands[i].setValue(v);
				intensities[i] = v;
				++i;
			}
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
}

void Spectrum::setup(Spectrum& spec) {
	// Instantiate the bands on the new spectrum.
	for(const Band& b : bands)
		spec.bands.emplace_back(b.wl(), 0);
}

void Spectrum::reset() {
	for(Band& b : bands) {
		b.setValue(0);
		b.reset();
	}
}

const std::string& Spectrum::projection() const {
	return m_projection;
}

void Spectrum::transform(double* trans) const {
	for(int i = 0; i < 6; ++i)
		trans[i] = m_trans[i];
}

std::unique_ptr<GDALReader>& Spectrum::raster() {
	return m_raster;
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

void Spectrum::write(GDALWriter* wtr, double minWl, double maxWl, int col, int row) {
	std::vector<double> v;
	for(const Band& b : bands) {
		if(b.wl() >= minWl && b.wl() <= maxWl)
			v.push_back(b.value());
	}
	wtr->write(v, col, row, 1, 1);
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
		m_bandProps.emplace(std::piecewise_construct, std::forward_as_tuple(band), std::forward_as_tuple(band, wl, fwhm));
		if(wl < minWl) minWl = wl;
		if(wl > maxWl) maxWl = wl;
	}
}

const std::map<int, BandProp>& BandPropsReader::bands() const {
	return m_bandProps;
}

void BandPropsReader::configureSpectrum(Spectrum& spec) {
	// Resize the spectrum to the number of bands in the configuration.
	spec.bands.resize(m_bandProps.size());
	spec.wavelengths.resize(spec.bands.size());
	spec.intensities.resize(spec.bands.size());
	// Set all the wavelengths on the bands.
	size_t i = 0;
	for(const auto& p : m_bandProps) {
		spec.bands[i].setWl(p.second.wl);
		spec.wavelengths[i] = p.second.wl;
		++i;
	}
}

class BinTree {
private:
	BinTree* left;
	BinTree* right;
	const Kernel* kernel;
	const Kernel& closer(double key, const Kernel& a, const Kernel& b) {
		if(std::abs(key - a.wl()) < std::abs(key - b.wl())) {
			return a;
		} else {
			return b;
		}
	}

	void init(const std::vector<Kernel>& k, int start, int end) {
		int mid = (end + start) / 2;
		kernel = &k[mid];
		if(mid > start) {
			left = new BinTree();
			left->init(k, start, mid);
		}
		if(end > mid + 1) {
			right = new BinTree();
			right->init(k, mid + 1, end);
		}
	}

public:

	BinTree() :
		left(nullptr), right(nullptr), kernel(nullptr) {
	}

	void init(const std::vector<Kernel>& k) {
		init(k, 0, (int) k.size());
	}

	const Kernel& find(double key) {
		if(key == kernel->wl()) {
			return *kernel;
		} else if(key < kernel->wl()) {
			return left ? closer(key, *kernel, left->find(key)) : *kernel;
		} else {
			return right ? closer(key, *kernel, right->find(key)) : *kernel;
		}
	}

	~BinTree() {
		if(left) delete left;
		if(right) delete right;
	}


};

void Convolve::run(ConvolveListener& listener,
		const std::string& bandDef, const std::string& bandDefDelim,
		const std::string& spectra, const std::string& spectraDelim,
		int spectraFirstRow, int spectraFirstCol,
		int spectraDateCol, int spectraTimeCol,
		const std::string& output, const std::string& outputDelim, FileType outputType,
		double inputScale, double tolerance, double bandShift, size_t memLimit, bool& running) {

	// Notify a listener.
	listener.started(this);

	// Load the band properties.
	BandPropsReader rdr;
	rdr.load(bandDef, bandDefDelim);

	const std::map<int, BandProp>& bands = rdr.bands();

	// Load the spectrum.
	Spectrum spec(spectraFirstRow, spectraFirstCol, spectraDateCol, spectraTimeCol);
	spec.load(spectra, spectraDelim, memLimit);
	spec.shift(bandShift);
	spec.scale(inputScale);

	// Configure the output (convolved) spectrum.
	Spectrum out;
	rdr.configureSpectrum(out);

	// Configure the kernels.
	int windowSize = 15;
	BinTree tree;
	std::vector<Kernel> kernels;
	{
		int i = 0;
		for(const auto& it : bands)
			kernels.emplace_back(it.second.wl, it.second.fwhm, windowSize, i++);
		std::sort(kernels.begin(), kernels.end());
		tree.init(kernels);
	}

	FileType ftype = getFileType(output);
	std::unique_ptr<Writer> writer;
	if(ftype == FileType::CSV) {
		writer.reset(new CSVWriter(output)); // wavelengths, bandNames
	} else {
		writer.reset(new GDALWriter(output, FileType::ENVI, spec.raster()->cols(), spec.raster()->rows(), rdr.bands().size())); //, wavelengths, bandNames
		static_cast<GDALWriter*>(writer.get())->setProjection(spec.projection());
		double trans[6];
		spec.transform(trans);
		static_cast<GDALWriter*>(writer.get())->setTransform(trans);
	}

	// Run the convolution record-by-record.
	size_t complete = 0, count = spec.count(); // Status counters.
	bool header = false;
	char delim = outputDelim[0];

	while(running && spec.next()) {
		out.reset();
		out.date = spec.date;
		out.time = spec.time;
		for(size_t i = 0; i < spec.wavelengths.size(); ++i) {
			const Kernel& k = tree.find(spec.wavelengths[i]);
			out.intensities[k.index()] = k.apply(spec.intensities, spec.wavelengths, i);
			if(!running) break;
		}
		out.update();

		if(!header && ftype == FileType::CSV) {
			// Write the header if this is the first iteration.
			out.writeHeader(static_cast<CSVWriter*>(writer.get())->outstr(), rdr.minWl, rdr.maxWl, delim);
			header = true;
		}

		// Write the record.
		if(ftype == FileType::CSV) {
			out.write(static_cast<CSVWriter*>(writer.get())->outstr(), rdr.minWl, rdr.maxWl, delim);
		} else {
			out.write(static_cast<GDALWriter*>(writer.get()), rdr.minWl, rdr.maxWl, spec.col(), spec.row());
		}

		// Update progress.
		m_progress = (double) complete++ / count;
		listener.update(this);
		if(!running) break;
	}
	// Notify listener of completion.
	listener.finished(this);
}

size_t Convolve::guess(const std::string& spectra, const std::string& spectraDelim) {
	Spectrum spec;
	return spec.inputSize(spectra, spectraDelim) * 1.5;
}

double Convolve::progress() const {
	return m_progress;
}

