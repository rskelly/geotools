//============================================================================
// Name        : contrem.cpp
// Author      : Rob Skelly
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <unistd.h>
#include <algorithm>

#include <gdal_priv.h>

#include <QtCore/QObject>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>

#include "ui/contrem_ui.hpp"
#include "contrem.hpp"
#include "reader.hpp"
#include "writer.hpp"
#include "plotter.hpp"

using namespace hlrg;
using namespace geo::plot;

class ContremApp : public QApplication {
public:
	ContremApp(int &argc, char **argv) :
		QApplication(argc, argv) {}

	bool notify(QObject *receiver, QEvent *e) {
		try {
			//if(plotter)
			//	plotter->step();
			return QApplication::notify(receiver, e);
		} catch(const std::exception &ex) {
			QMessageBox err;
			err.setText("Error");
			err.setInformativeText(QString(ex.what()));
			err.exec();
			return false;
		}
	}
};

int runWithGui(ContremApp& app, ContremForm& form) {
	QDialog qform;
	form.setupUi(&qform);
	qform.show();
	return app.exec();
}

class DummyListener : public ContremListener {
public:
	void started(Contrem*) {
		std::cout << "Started\n";
	}
	void update(Contrem* conv) {
		std::cout << "Progress: " << (conv->progress() * 100) << "%\n";
	}
	void stopped(Contrem*) {
		std::cout << "Stopped.\n";
	}
	void finished(Contrem*) {
		std::cout << "Finished.\n";
	}
};

void usage() {
	std::cerr << "Usage: contrem [options]\n"
			<< " -sf A GDAL-readable data file containing spectral samples; can contain any number of bands >= 2.\n"
			<< " -st The spectrum file type; GTiff, ENVI or CSV.\n"
			<< " -rf A GDAL-readable mask file.\n"
			<< " -bf A CSV file containing a mapping from wavelength to (1-based) band index.\n"
			<< " -of An output file template. This is a filename with no extension that will be modified as\n"
			<< "     appropriate. Parent directories will be created.\n"
			<< " -od The driver to use for output rasters. GTiff or ENVI.\n"
			<< " -oe File extension for raster files. Defaults to .dat for ENVI files, .tif for GTiff.\n"
			<< " -w  An integer giving the (0-based) column index in -b which contains wavelengths.\n"
			<< " -i  An integer giving the (0-based) column index in -b which contains the band indices.\n"
			<< " -z  If given, indicates the presence of a header in the band map that must be skipped.\n"
			<< " -l  The minimum wavelength to consider.\n"
			<< " -h  The maximum wavelength to consider.\n"
			<< " -t  The number of threads to use. Default 2.\n"
			<< " -nm Normalization method. ConvexHull, ConvexHullLongestSeg or Line.\n"
			<< "Run without arguments for GUI version.\n";
}

int main(int argc, char** argv) {

	int ret = 0;

	Plotter plot; // Will start if any jobs are sent to it.

	if(argc > 1) {

		try {
			Contrem contrem;

			contrem.extension = ".dat";
			contrem.threads = 2;
			contrem.spectraType = FileType::Unknown;
			contrem.normMethod = NormMethod::Unknown;
			contrem.outputType = FileType::Unknown;

			for(int i = 0; i < argc; ++i) {
				std::string arg(argv[i]);
				if(arg == "-sf") {
					contrem.spectra = argv[++i];
				} else if(arg == "-st") {
					std::string s(argv[++i]);
					if(s == "CSV") {
						contrem.spectraType = FileType::CSV;
					} else if(s == "GTiff") {
						contrem.spectraType = FileType::GTiff;
					} else if(s == "ENVI") {
						contrem.spectraType = FileType::ENVI;
					}
				} else if(arg == "-rf") {
					contrem.roi = argv[++i];
				} else if (arg == "-bf") {
					// contrem.bandFile = argv[++i];
					std::cerr << "Band file not implemented.\n";
					return 1;
				} else if(arg == "-l") {
					contrem.minWl = atof(argv[++i]);
				} else if(arg == "-h") {
					contrem.maxWl = atof(argv[++i]);
				} else if(arg == "-w") {
					//wlCol = atoi(argv[++i]);
				} else if(arg == "-i") {
					//bandCol = atoi(argv[++i]);
				} else if(arg == "-z") {
					//bandHeader = true;
				} else if(arg == "-of") {
					contrem.output = argv[++i];
				} else if(arg == "-oe") {
					contrem.extension = argv[++i];
				} else if(arg == "-od") {
					std::string d(argv[++i]);
					if(d == "GTiff") {
						contrem.outputType = FileType::GTiff;
					} else if(d == "ENVI") {
						contrem.outputType = FileType::ENVI;
					} else if(d == "CSV") {
						contrem.outputType = FileType::CSV;
					}
				} else if(arg == "-t") {
					contrem.threads = atoi(argv[++i]);
				} else if(arg == "-nm") {
					std::string d(argv[++i]);
					if(d == "ConvexHull") {
						contrem.normMethod = NormMethod::ConvexHull;
					} else if(d == "ConvexHulLongestSeg") {
						contrem.normMethod = NormMethod::ConvexHullLongestSeg;
					} else if(d == "Line") {
						contrem.normMethod = NormMethod::Line;
					}
				}
			}

			DummyListener dl;

			contrem.running = true;
			contrem.run(&dl);

		} catch(const std::exception& ex) {
			std::cerr << ex.what() << "\n";
			usage();
			ret = 0;
		}

	} else {

		ContremApp app(argc, argv);
		ContremForm form(&app);

		runWithGui(app, form);

	}

	plot.stop();

	return ret;
}
