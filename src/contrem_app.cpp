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

int main(int argc, char** argv) {

	Plotter plot;
	if(!plot.start())
		return 1;

	ContremApp app(argc, argv);
	ContremForm form(&app);

	runWithGui(app, form);

	plot.stop();

	return 0;
}


/*



class DummyListener : public ContremListener {
public:
	void started(Contrem* conv) {
		std::cout << "Started\n";
	}
	void update(Contrem* conv) {
		std::cout << "Progress: " << (conv->progress() * 100) << "%\n";
	}
	void stopped(Contrem* conv) {
		std::cout << "Stopped.\n";
	}
	void finished(Contrem* conv) {
		std::cout << "Finished.\n";
	}
};

void usage() {
	std::cerr << "Usage: contrem [options]\n"
			<< " -d A GDAL-readable data file containing spectral samples; can contain any number of bands >= 2.\n"
			<< " -r An ENVI ROI text file.\n"
			<< " -b A CSV file containing a mapping from wavelength to (1-based) band index.\n"
			<< " -w An integer giving the (0-based) column index in -b which contains wavelengths.\n"
			<< " -i An integer giving the (0-based) column index in -b which contains the band indices.\n"
			<< " -z If given, indicates the presence of a header in the band map that must be skipped.\n"
			<< " -o An output file template. This is a filename with no extension that will be modified as\n"
			<< "    appropriate. Parent directories will be created.\n"
			<< " -l The minimum wavelength to consider.\n"
			<< " -h The maximum wavelength to consider.\n"
			<< " -s The size of the buffer. Default is 256. Larger buffers are possible, but one must\n"
			<< "    consider that multiple buffers may be in memory at once.\n"
			<< " -t The number of threads to use. Default 2.\n"
			<< " -p By default, sample statistics are used. This flag forces the use of\n"
			<< "    population statistics.\n"
			<< " -v The driver to use for output rasters. Defaults to ENVI, but any GDAL-writable\n"
			<< "    format will do.\n"
			<< " -e File extension for raster files. Defaults to .dat for ENVI files.\n";
}

int main(int argc, char** argv) {
*/
	/*
	if(argc > 1) {

		Contrem processor;

		processor.bufferSize = 256;
		processor.sampleStats = false;
		processor.driver = "ENVI";
		processor.extension = ".dat";
		processor.threads = 1;

		int wlCol = -1;
		int bandCol = -1;
		bool bandHeader = false;
		double minWl = 0;
		double maxWl = 0;
		std::string datafile;
		std::string roifile;
		std::string bandfile;

		try {
			int c;
			while((c = getopt(argc, argv, "d:r:b:o:l:h:s:w:i:t:e:v:zp")) != -1) {
				switch(c) {
				case 'd': datafile = optarg; break;
				case 'r': roifile = optarg; break;
				case 'b': bandfile = optarg; break;
				case 'l': minWl = atof(optarg); break;
				case 'h': maxWl = atof(optarg); break;
				case 'w': wlCol = atoi(optarg); break;
				case 'i': bandCol = atoi(optarg); break;
				case 'z': bandHeader = true; break;
				case 'o': processor.outfile = optarg; break;
				case 's': processor.bufferSize = atoi(optarg); break;
				case 't': processor.threads = atoi(optarg); break;
				case 'p': processor.sampleStats = false; break;
				case 'v': processor.driver = optarg; break;
				case 'e': processor.extension = optarg; break;
				default: break;
				}
			}

			if(processor.bufferSize <= 0)
				throw std::invalid_argument("Buffer size must be larger than zero.");
			if(datafile.empty() && roifile.empty())
				throw std::invalid_argument("Data  or ROI file not given.");
			if(processor.outfile.empty())
				throw std::invalid_argument("Output file template not given.");
			if(processor.threads < 1)
				throw std::invalid_argument("At least one thread is required.");
			if(!bandfile.empty() && (bandCol < 0 || wlCol < 0))
				throw std::invalid_argument("If the band file is given, wavelength and band columns must also be given.");

			Reader* reader;
			if(!roifile.empty()) {
				reader = new ROIReader(roifile);
			} else if(!datafile.empty()) {
				reader = new GDALReader(datafile);
			} else {
				throw std::invalid_argument("No input file (-r or -d) given.");
			}

			if(!bandfile.empty()) {
				if(wlCol == -1 || bandCol == -1 || wlCol == bandCol)
					throw std::invalid_argument("If a band file is given, the indices for "
							"wavelength and band must be >=0 and different from each other.");
				BandMapReader br(bandfile, wlCol, bandCol, bandHeader);
				reader->setBandMap(br.bandMap());
			}

			if(minWl > 0 && maxWl > 0)
				reader->setBandRange(minWl, maxWl);

			reader->setBufSize(processor.bufferSize);

			DummyListener dl;

			processor.run(&dl, reader);

		} catch(const std::exception& ex) {
			std::cerr << ex.what() << "\n";
			usage();
			return 1;
		}

	} else {
*/
