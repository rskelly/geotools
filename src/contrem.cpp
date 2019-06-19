/*
 * processor.cpp
 *
 *  Created on: May 9, 2018
 *      Author: rob
 */

#include <vector>
#include <list>
#include <iostream>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <thread>
#include <condition_variable>
#include <algorithm>

#include <geos_c.h>

#include "util.hpp"
#include "contrem.hpp"
#include "reader.hpp"
#include "writer.hpp"

using namespace hlrg::contrem;
using namespace hlrg::reader;
using namespace hlrg::writer;

namespace {

	std::unique_ptr<Reader> getReader(const std::string& file, bool transpose, int headerRows, int minCol, int maxCol, int idCol) {
		std::unique_ptr<Reader> rdr;
		FileType type = getFileType(file);
		switch(type) {
		case FileType::CSV:
			rdr.reset(new CSVReader(file, transpose, headerRows, minCol, maxCol, idCol));
			break;
		case FileType::GTiff:
		case FileType::ENVI:
			rdr.reset(new GDALReader(file));
			break;
		case FileType::ROI:
		case FileType::SHP:
		case FileType::SQLITE:
		case FileType::Unknown:
		default:
			throw std::runtime_error("Unknown file type for " + file);
		}
		return rdr;
	}

	/**
	 * A line, representing a segment from a convex hull.
	 */
	class line {
	public:
		double x0;
		double y0;
		double x1;
		double y1;

		line() : line(0, 0, 0, 0) {}

		line(double x0, double y0, double x1, double y1) :
			x0(x0), y0(y0), x1(x1), y1(y1) {}

		double length() const {
			return std::sqrt(std::pow(x0 - x1, 2.0) + std::pow(y0 - y1, 2.0));
		}

		double slope() const {
			return (y1 - y0) / (x1 - x0);
		}

		double yint() const {
			return -slope() * x0;
		}
	};

	/**
	 * An input point. Contains a wavelength and sample intensity.
	 */
	class inpoint {
	public:
		double w; // Wavelength
		double ss; // Sample spectra
		inpoint(double w, double ss) :
			w(w), ss(ss) {}
	};

	/**
	 * An output point. Contains the same info as an input point,
	 * plus continuum removal values.
	 */
	class outpoint {
	public:
		double w; 		// Wavelength.
		double ss; 		// Sample spectra (intensity).
		double ch; 		// Intersection with convex hull (y).
		double cr; 		// Continuum removal (ss/ch).
		double crn;		// Continuum removal normalized (1 - cr).
		double crm;		// Mirrored cr.
		double crnm;	// Mirrored normalized cr.

		outpoint(inpoint& in, double ch) :
			w(in.w), ss(in.ss),
			ch(ch), cr(0), crn(0), crm(0), crnm(0) {}
		outpoint(inpoint& in) :
			outpoint(in, 0) {}
	};

	/**
	 * An input pixel, which includes spectral information
	 * and cell coordinates.
	 */
	class input {
	public:
		std::string id;		///<! Identifies the datum. Useful for spreadsheets.
		int c, r;
		std::vector<inpoint> data;

		input() : input("", 0, 0) {}

		input(const std::string& id, int c, int r) :
			id(id),
			c(c), r(r) {
		}

	};

	/**
	 * An output pixel which contains convex hull information,
	 * continuum removal information and cell coordinates.
	 */
	class output {
	public:
		std::string id;					///<! Identifies the datum. Useful for spreadsheets.
		int c, r;						///<! The column and row of the datum. Useful for rasters.
		double area; 					///<! Total hull area.
		double larea; 					///<! Left hull area.
		double rarea; 					///<! Right hull area.
		double symmetry; 				///<! Ratio of left/right hull areas.
		double maxDepth; 				///<! Maximum depth w/r/t line.
		double maxWl; 					///<! Wavelength at the maxCrm value.
		double slope;					///<! The slope of the regression line through the crnm points.
		double yint;					///<! The y-intercept of the regression.
		int maxCount; 					///<! The number of equal maximum values.
		int maxIdx;						///<! The index of the first maximum.

		std::vector<outpoint> data;		///<! The output data.
		std::vector<double> hullx;
		std::vector<double> hully;

		output() : output("", 0, 0) {}

		output(input& in) :
			output(in.id, in.c, in.r) {}

		output(const std::string& id, int c, int r) :
			id(id),
			c(c), r(r),
			area(0), larea(0), rarea(0), symmetry(0),
			maxDepth(0), maxWl(0),
			slope(0), yint(0),
			maxCount(0),
			maxIdx(0) {}

		bool compute(const std::vector<line>& lines) {
			maxIdx = 0;
			maxDepth = 0;
			int i = 0;
			double dif;
			for(outpoint& pt : data) {
				pt.cr = pt.ss / pt.ch;		// Continuum removal -- intensity proportional to height of hull.
				pt.crm = 1 - pt.cr;			// Mirrored continuum removal.
				if((dif = pt.ch - pt.ss) > maxDepth) {
					maxDepth = dif;
					maxIdx = i;
					maxWl = pt.w;
					++maxCount;
				} else if(dif == maxDepth) {
					++maxCount;
				}
				++i;
			}

			if(maxDepth <= 0)
				return false;

			double sxy = 0, sx = 0, sy = 0, sxx = 0;
			size_t n = data.size();

			area = 0;
			larea = 0;

			// Calculate the other ch metrics.
			for(size_t i = 0; i < n; ++i) {
				outpoint& pt = data[i];
				pt.crn = pt.crm / maxDepth;		// Normalized continuum removal -- normalized against maximum depth.
				pt.crnm = 1 - pt.crn;			// Mirrored normalized continuum removal.
				if(i > 0 && i < n - 1) {
					sxy += pt.w * pt.crnm;
					sx += pt.w;
					sy += pt.crnm;
					sxx += pt.w * pt.w;
				}

				// Compute the left and overall spectrum area.
				if(i < n - 1) {
					double ww = data[i + 1].w - data[i].w;
					double a = ((data[i].crn + data[i + 1].crn) * ww) / 2.0; // Add the heights, times the width, divide by two.
					area += a;
					if((int) i <= maxIdx)
						larea += a;
				}
			}

			// Compute the rest of the numbers.
			if(area == 0 || larea == 0 || larea == area) {
				// If the left or right area is zero, or the total is zero, we have a problem.
				larea = 0;
				area = 0;
				symmetry = 0;
				rarea = 0;
			} else {
				rarea = area - larea;
				symmetry = larea / rarea;
			}

			double maxLen = 0, len;
			size_t lineIdx = 0;
			for(size_t i = 0; i < lines.size(); ++i) {
				if((len = lines[i].length()) > maxLen) {
					maxLen = len;
					lineIdx = i;
				}
			}

			slope = lines[lineIdx].slope();
			yint = lines[lineIdx].yint();

			return true;
		}
	};

	/**
	 * Returns the y value corresponding to x along the given line.
	 * NaN if no intersection found.
	 */
	double interpolate(double x, double x0, double y0, double x1, double y1) {
		if(x < x0 || x > x1) {
			return std::numeric_limits<double>::quiet_NaN();
		} else if(x == x0) {
			return y0;
		} else if(x == x1) {
			return y1;
		} else {
			return y0 + (x - x0) / (x1 - x0) * (y1 - y0);
		}
	}

	/**
	 * Compute the convex hull around the points and return the line segments.
	 */
	std::vector<line> convexHull(const std::vector<inpoint>& in, double* area = nullptr) {

		// Make a list of Coordinates.
		GEOSCoordSequence* seq;
		std::vector<GEOSGeometry*> coords;
		for(const inpoint& pt : in) {
			seq = GEOSCoordSeq_create(1, 2);
			GEOSCoordSeq_setX(seq, 0, pt.w);
			GEOSCoordSeq_setY(seq, 0, pt.ss);
			coords.push_back(GEOSGeom_createPoint(seq));
		}

		// Make a MultiPoint from the coords and get the ConvexHull.
		GEOSGeometry* mp = GEOSGeom_createCollection(GEOS_MULTIPOINT, coords.data(), coords.size());
		GEOSGeometry* hull = GEOSConvexHull(mp);

		// Get the area for output.
		if(area != nullptr)
			*area = GEOSArea(hull, area);

		// Extract the line segments.
		std::vector<line> lines;
		const GEOSGeometry* ring = GEOSGetExteriorRing(hull);
		GEOSGeometry* p0, *p1;
		double x0, x1, y0, y1;
		for(int i = 0; i < GEOSGeomGetNumPoints(ring) - 1; ++i) {
			p0 = GEOSGeomGetPointN(ring, i);
			p1 = GEOSGeomGetPointN(ring, i + 1);
			GEOSGeomGetX(p0, &x0);
			GEOSGeomGetY(p0, &y0);
			GEOSGeomGetX(p1, &x1);
			GEOSGeomGetY(p1, &y1);
			// Only add a segment if it isn't a bottom or end segment (which have at least one zero y-coordinate).
			if(y0 > 0 && y1 > 0)
				lines.emplace_back(x0, y0, x1, y1);

			GEOSGeom_destroy(p0);
			GEOSGeom_destroy(p1);
		}

		GEOSGeom_destroy(mp);
		GEOSGeom_destroy(hull);

		return lines;
	}

	class QConfig {
	private:
		int steps;							///<! The total number of steps to complete the processing. Used for the status bar.
		int step;							///<! The current completed step.

	public:
		Plotter* plotter;
		Contrem* contrem;
		std::list<input> inqueue;
		std::list<output> outqueue;
		std::mutex inmtx;
		std::mutex outmtx;
		std::mutex readmtx;
		std::condition_variable incv;
		std::condition_variable outcv;
		std::condition_variable readcv;
		bool inRunning;
		bool outRunning;
		bool useROI;						///<! If the spectra file is the right type, the ROI can be used.
		int cols;
		int rows;
		int bands;

		std::vector<double> wavelengths;
		std::vector<std::string> bandNames;

		std::vector<std::string> getWavelengthNames() const {
			std::vector<std::string> names;
			for(double w : wavelengths)
				names.push_back(std::to_string(w));
			return names;
		}
	};


	/**
	 * Return a list of lines depending on configuration.
	 * 1) If a hull is desired.
	 * 2) If the longest segment of a hull is desired.
	 * 3) If a straight line from start to end of spectrum.
	 *
	 * \param config The config object.
	 * \param pts The list of input points.
	 * \return A list of lines.
	 */
	std::vector<line> getLines(QConfig* config, std::vector<inpoint>& pts) {

		std::vector<line> lines;

		NormMethod method = config->contrem->normMethod;
		if(method == NormMethod::ConvexHull || method == NormMethod::ConvexHullLongestSeg) {

			std::vector<inpoint> _pts(pts);
			_pts.emplace_back(pts.back().w, 0.0);
			_pts.emplace_back(pts.front().w, 0.0);

			// Compute the hull.
			lines = convexHull(_pts);

			if(method == NormMethod::ConvexHullLongestSeg) {
				// If required, take the longest segment out of the hull and use it for normalization.
				size_t idx = 0;
				double len = 0;
				for(size_t i = 0; i < lines.size(); ++i) {
					if(lines[i].length() > len) {
						len = lines[i].length();
						idx = i;
					}
				}
				lines[0] = std::move(lines[idx]);
				lines.resize(1);
			}

		} else {

			// Use a single line from the first to the last point on the spectrum.
			lines.emplace_back(pts.front().w, pts.front().ss, pts.back().w, pts.back().ss);

		}

		return lines;
	}


	/**
	 * Process the input queue.
	 */
	void processQueue(QConfig* config) {

		initGEOS(0, 0);

		input in;
		while(config->contrem->running) {
			std::this_thread::yield();
			{
				std::unique_lock<std::mutex> lk(config->inmtx);
				// If the input is empty or the output is too large, pause.
				while((config->inqueue.empty() || config->outqueue.size() > 1000) && config->inRunning)
					config->incv.wait(lk);
				// If input is empty and reading is done, quit the loop.
				if(config->inqueue.empty() && !config->inRunning)
					break;
				in = config->inqueue.front();
				config->inqueue.pop_front();
			}

			if(!config->contrem->running)
				break;

			// Adjust <=0 intensities to MIN_VALUE. This enables the creation
			// of a hull even though the area of the hull will be zero for practical purposes.
			for(size_t i = 0; i < in.data.size(); ++i) {
				if(in.data[i].ss <= MIN_VALUE)
					in.data[i].ss = MIN_VALUE;
			}

			// Add two corner points to complete the hull.
			std::vector<inpoint> pts(in.data.begin(), in.data.end());

			// Get the linework for normalization.
			std::vector<line> lines = getLines(config, pts);

			// Create an output pixel to hold computed values.
			output out(in);


			if(!config->contrem->running)
				break;

			// Store the points of the convex hull/line in case plotting is desired.
			if(config->contrem->plotOrig) {
				out.hullx.clear();
				out.hully.clear();
				for(line& l : lines) {
					out.hullx.push_back(l.x0);
					out.hully.push_back(l.y0);
				}
				out.hullx.push_back(lines.back().x1);
				out.hully.push_back(lines.back().y1);
			}

			// Find the intersection point for each wavelength.
			// If an intersection isn't found, discard the point (this may
			// occur if the single line from the hull is used.
			for(inpoint& pt : pts) {
				for(line& l : lines) {
					double ch = interpolate(pt.w, l.x0, l.y0, l.x1, l.y1);
					if(!std::isnan(ch) && pt.ss != 0) {
						out.data.emplace_back(pt, ch);
						break;
					}
				}
			}

			if(out.data.size() < 2)
				throw std::runtime_error("The list of input points is too small.");

			if(!config->contrem->running)
				break;

			// We were going to do interpolation for adjacent maxima, but put it off.
			// Kopăcková, V., & Koucká, L. (2017). Integration of absorption feature information from visible to
			// longwave infrared spectral ranges for mineral mapping. Remote Sensing, 9(10), 8–13. https://doi.org/10.3390/rs9101006
			// If there are  >2 maxima, or the distance between them is > than the configured
			// interp distance, flag the cell and move on. Otherwise, interpolate.

			// Calculate the cr and crm, etc., and get the max value and index.
			if(out.compute(lines)) {

				// Send to output queue and notify.
				std::lock_guard<std::mutex> lk(config->outmtx);
				config->outqueue.push_back(out);
			}

			config->contrem->nextStep();

			config->outcv.notify_one();
			config->readcv.notify_one();

		}

		finishGEOS();

	}



	/**
	 * Process the output queue and write to file.
	 */
	void writeQueue(QConfig* config) {

		std::string outfile = config->contrem->output;
		FileType outfileType = config->contrem->outputType;
		const std::vector<double>& wavelengths = config->wavelengths;
		const std::vector<std::string>& bandNames = config->bandNames;
		int cols = config->cols;
		int rows = config->rows;
		int bands = config->bands;

		std::string ext;
		switch(outfileType) {
		case FileType::ENVI:
			ext = "";
			break;
		case FileType::GTiff:
			ext = ".tif";
			break;
		case FileType::CSV:
			ext = ".csv";
			break;
		default:
			throw std::invalid_argument("Unknown output type: " + fileTypeAsString(outfileType));
		}

		// Remove the extension if there is one.
		outfile = outfile.substr(0, outfile.find_last_of("."));
		std::string outdir = outfile.substr(0, outfile.find_last_of("/"));
		std::string plotdir = outdir + "/hull_img";

		if(isfile(outdir))
			throw std::runtime_error("The output directory is an extant file.");
		if(!isdir(outdir) && !makedir(outdir))
			throw std::runtime_error("Failed to create output directory.");
		if(!isdir(plotdir) && !makedir(plotdir))
			throw std::runtime_error("Failed to create plot directory.");

		if(!config->contrem->running)
			return;

		char* meta;
		std::unordered_map<std::string, std::unique_ptr<Writer>> writer;

		if(outfileType == FileType::CSV) {
			writer["writerss"].reset(new CSVWriter(outfile + "_ss" + ext, wavelengths, bandNames));
			writer["writerch"].reset(new CSVWriter(outfile + "_ch" + ext, wavelengths, bandNames));
			writer["writercr"].reset(new CSVWriter(outfile + "_cr" + ext, wavelengths, bandNames));
			writer["writercrnm"].reset(new CSVWriter(outfile + "_crnm" + ext, wavelengths, bandNames));
			writer["writerhull"].reset(new CSVWriter(outfile + "_agg" + ext, {}, {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crm", "max_crm_wl", "max_count", "slope", "y-int"}));
			writer["writermax"].reset(new CSVWriter(outfile + "_maxcount" + ext, {}, {"equal_max_count"}));
			writer["writervalid"].reset(new CSVWriter(outfile + "_valid" + ext, {}, {"valid_hull"}));
		} else {
			writer["writerss"].reset(new GDALWriter(outfile + "_ss" + ext, outfileType, cols, rows, bands, wavelengths, bandNames));
			writer["writerch"].reset(new GDALWriter(outfile + "_ch" + ext, outfileType, cols, rows, bands, wavelengths, bandNames));
			writer["writercr"].reset(new GDALWriter(outfile + "_cr" + ext, outfileType, cols, rows, bands, wavelengths, bandNames));
			writer["writercrnm"].reset(new GDALWriter(outfile + "_crnm" + ext, outfileType, cols, rows, bands, wavelengths, bandNames));
			writer["writerhull"].reset(new GDALWriter(outfile + "_agg" + ext, outfileType, cols, rows, 9, {}, {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crm", "max_crm_wl", "max_count", "slope", "y-int"}));
			writer["writermax"].reset(new GDALWriter(outfile + "_maxcount" + ext, outfileType, cols, rows, 1, {}, {"equal_max_count"}, &meta, DataType::Byte));
			writer["writervalid"].reset(new GDALWriter(outfile + "_valid" + ext, outfileType, cols, rows, 1, {}, {"valid_hull"}, &meta, DataType::Byte));
			writer["writermax"]->fill(0);
			writer["writervalid"]->fill(0);
		}

		// Initialize the output lists.
		std::vector<double> ss;
		std::vector<double> ch;
		std::vector<double> cr;
		std::vector<double> crnm;
		std::vector<double> hull;
		std::vector<double> w;
		std::vector<int> maxima(1);
		std::vector<int> valid(1);

		std::fill(maxima.begin(), maxima.end(), 0);
		std::fill(valid.begin(), valid.end(), 0);

		// If there's a sample points file and a raster reader, we can use the sample points.
		std::unique_ptr<PointSetReader> samples;
		bool hasSamples = false;
		if(!config->contrem->samplePoints.empty() && config->contrem->grdr) {
			samples.reset(new PointSetReader(config->contrem->samplePoints,
					config->contrem->samplePointsLayer, config->contrem->samplePointsIDField));
			samples->toGridSpace(config->contrem->grdr);
			hasSamples = true;
		}

		// Keeps track of the number of unique completed items for progress tracking.
		int count = 0;

		// Processor loop.
		output out;
		while(config->contrem->running) {

			{
				// Get an item from the queue.
				std::unique_lock<std::mutex> lk(config->outmtx);
				while(config->outqueue.empty() && config->outRunning)
					config->outcv.wait(lk);
				if(config->outqueue.empty() && !config->outRunning)
					break;
				out = config->outqueue.front();
				config->outqueue.pop_front();
			}

			++count;

			if(!config->contrem->running)
				break;

			// Populate the lists for plotting/aggregating.
			for(const outpoint& o : out.data) {
				ss.push_back(o.ss);
				ch.push_back(o.ch);
				cr.push_back(o.cr);
				crnm.push_back(o.crnm);
				w.push_back(o.w);
			}

			if(!config->contrem->running)
				break;

			if(hasSamples) {
				hlrg::reader::Point pt;
				pt.c(out.c);
				pt.r(out.r);
				if(samples->sampleNear(pt, 0.5)) {

					// If appropriate plot the normalized spectrum.
					if(config->contrem->plotNorm){
						std::string plotfile = plotdir + "/norm_" + sanitize(out.id) + "_" + std::to_string(out.c) + "_" + std::to_string(out.r) + ".png";
						std::string title = "Normalized Spectrum (" + pt.id() + ", " + std::to_string(pt.x()) + "," + std::to_string(pt.x()) + ")";
						std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> items;
						items.emplace_back("Normalized Spectrum", w, crnm);
						//items.emplace_back("Regression", std::vector<double>({w.front(), w.back()}), std::vector<double>({w.front() * out.slope + out.yint, w.back() * out.slope + out.yint}));
						Plotter::instance().enqueue(plotfile, title, items);
					}
					if(config->contrem->plotOrig){
						std::string plotfile = plotdir + "/orig_" + sanitize(out.id) + "_" + std::to_string(out.c) + "_" + std::to_string(out.r) + ".png";
						std::string title = "Original Spectrum + Hull (" + pt.id() + ", " + std::to_string(pt.x()) + "," + std::to_string(pt.y()) + ")";
						std::vector<std::tuple<std::string, std::vector<double>, std::vector<double>>> items;
						items.emplace_back("Original Spectrum", w, ss);
						items.emplace_back("Convex Hull", out.hullx, out.hully);
						Plotter::instance().enqueue(plotfile, title, items);
					}
				}
			}

			if(!config->contrem->running)
				break;

			// The number of equal maxima.
			maxima[0] = out.maxCount > 1 ? 0 : 1;

			// A hull is valid if the area, left area and right area are non-zero.
			valid[0] = out.area > 0 && out.rarea > 0 && out.larea > 0;

			// Data for a single hull.
			hull = {out.area, out.larea, out.rarea, out.symmetry, out.maxDepth, out.maxWl, (double) out.maxCount, out.slope, out.yint};

			writer["writerss"]->write(ss, out.c, out.r, 1, 1, 1, 1, out.id);
			writer["writerch"]->write(ch, out.c, out.r, 1, 1, 1, 1, out.id);
			writer["writercr"]->write(cr, out.c, out.r, 1, 1, 1, 1, out.id);
			writer["writercrnm"]->write(crnm, out.c, out.r, 1, 1, 1, 1, out.id);
			writer["writerhull"]->write(hull, out.c, out.r, 1, 1, 1, 1, out.id);
			writer["writermax"]->write(maxima, out.c, out.r, 1, 1, 1, 1, out.id);
			writer["writervalid"]->write(valid, out.c, out.r, 1, 1, 1, 1, out.id);

			ss.clear();
			ch.clear();
			cr.clear();
			crnm.clear();
			w.clear();

			config->contrem->nextStep();

			config->incv.notify_one();
		}

		writer["writerhull"]->writeStats(outfile + "_agg_stats.csv", {"hull_area", "hull_left_area", "hull_right_area", "hull_symmetry", "max_crm", "max_crm_wl", "max_count", "slope", "yint"});
	}

}

void Contrem::run(ContremListener* listener) {

	if(!listener)
		throw std::runtime_error("A listener is required.");

	QConfig config;
	m_listener = listener;
	m_listener->started(this);

	initSteps(1, 100);

	std::unique_ptr<Reader> reader = getReader(spectra, wlTranspose, wlHeaderRows, wlMinCol, wlMaxCol, wlIDCol);
	reader->setBandRange(minWl, maxWl);
	if(reader->fileType() != FileType::CSV) {
		grdr = static_cast<GDALReader*>(reader.get());
		grdr->remap(minWl, maxWl);
	}

	config.contrem = this;
	config.cols = reader->cols();
	config.rows = reader->rows();
	config.bands = reader->bands();
	config.useROI = getFileType(spectra) != FileType::CSV;

	// A list of wavelengths.
	config.wavelengths = reader->getWavelengths();
	config.bandNames = reader->getBandNames();

	config.inRunning = true;
	config.outRunning = true;

	nextStep();

	// A list of wavelengths as strings for labelling.
	std::vector<std::string> wavelengthMeta;
	for(double w : config.wavelengths)
		wavelengthMeta.push_back(std::to_string(w));

	// Start the processing threads.
	std::list<std::thread> t0;
	for(int i = 0; i < threads; ++i)
		t0.emplace_back(processQueue, &config);

	// Start the output thread.
	std::thread t1(writeQueue, &config);

	nextStep();

	// Check for a mask file.
	bool hasRoi = false;
	std::vector<bool> mask;
	{
		if(config.useROI && !roi.empty() && isfile(roi)) {
			try {
				GDALReader rdr(roi);
				int col, row, cols = 1;
				std::string id;
				std::vector<double> buf(rdr.cols());
				mask.resize(rdr.cols() * rdr.rows());
				while(running && rdr.next(buf, 1, cols, col, row)) {
					for(int i = 0; i < cols; ++i)
						mask[row * cols + i] = buf[i] > 0;
				}
				hasRoi = true;
			} catch(const std::exception& ex) {
				std::cerr << "Could not open mask: " << ex.what() << "\n";
			}
		}
	}

	nextStep();

	// Determine the number of steps for status-keeping.
	{
		int steps;
		switch(getFileType(spectra)) {
		case FileType::CSV:
			steps = reader->rows();
			break;
		case FileType::GTiff:
		case FileType::ENVI:
			if(hasRoi) {
				steps = std::count(mask.begin(), mask.end(), true);
			} else {
				steps = reader->rows() * reader->cols();
			}
			break;
		default:
			break;
		}

		// Each datum goes through three steps; 1) adding to the queue; 2) processing; 3) writing.
		// There are 3 steps at the end.
		initSteps(1, steps * 3 + 2);
	}


	// A buffer for input data. Stores a row from a raster, or a single "pixel" from a table.
	std::vector<double> buf(config.cols * reader->bands());

	// Read through the buffer and populate the input queue.
	int bands = reader->bands();
	int cols, col, row;
	std::string id;
	while(running && reader->next(id, buf, cols, col, row)) {

		{
			// Read out the input objects and add to the queue.
			std::lock_guard<std::mutex> lk(config.inmtx);

			// If there's a mask, check it. Skip if necessary.
			if(hasRoi && !mask[row * cols + col])
				continue;

			input in(id, col, row);
			for(int b = 0; b < bands; ++b) {
				double v = buf[b];
				double w = config.wavelengths[b];
				in.data.emplace_back(w, v);
			}
			config.inqueue.push_back(in);

			nextStep();
		}

		// Notify the input processor.
		config.incv.notify_all();

		// If input queue gets too large, wait before adding more.
		std::unique_lock<std::mutex> lk(config.readmtx);
		while(config.inqueue.size() > 1000)
			config.readcv.wait(lk);
	}

	nextStep();

	// Let the processor threads finish.
	config.inRunning = false;
	config.incv.notify_all();
	for(std::thread& t : t0)
		t.join();

	nextStep();

	// Let the output thread finish.
	config.outRunning = false;
	config.outcv.notify_all();
	t1.join();

	nextStep();

	listener->finished(this);
}

void Contrem::initSteps(int step, int steps) {
	m_step = step;
	m_steps = steps;
	m_listener->update(this);
}

void Contrem::nextStep() {
	m_step++;
	m_listener->update(this);
}

double Contrem::progress() const {
	return (double) m_step / m_steps;
}
