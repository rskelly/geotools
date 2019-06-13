/*
 * convolve.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: rob
 */

#include <iostream>

#include <QtCore/QObject>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>

#include "convolver.hpp"
#include "ui/convolve_ui.hpp"

using namespace hlrg;

int runWithGui(int argc, char **argv) {
	class ConvolveApp : public QApplication {
	public:
		ConvolveApp(int &argc, char **argv) : QApplication(argc, argv) {}
		bool notify(QObject *receiver, QEvent *e) {
			try {
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

	ConvolveApp q(argc, argv);
	Convolver conv;
	ConvolveForm form(&conv, &q);
	QDialog qform;
	form.setupUi(&qform);
	qform.show();
	return q.exec();
}

class DummyListener : public ConvolverListener {
public:
	void started(Convolver* conv) {
		std::cout << "Running ";
	}
	void update(Convolver* conv) {
		int p = (int) (conv->progress() * 100);
		if(p % 25 == 0)
			std::cout << " " << p << "% ";
		if(p % 10 == 0)
			std::cout << ".";
	}
	void stopped(Convolver* conv) {
		std::cout << " Stopped.\n";
	}
	void finished(Convolver* conv) {
		std::cout << " Done.\n";
	}
};

void usage() {
	std::cerr << "Usage: convolve [[options] <band definition file> <spectra file> <output file>]\n"
			<< " -t		The threshold -- the gaussian will extend out until it is below this value. (Default 0.0001).\n"
			<< " -f 	The wavelength shift. Used to counter mis-calibration in the spectrometer. (Default 0).\n"
			<< " -s		Scale the input. (Default 1).\n"
			<< " -ds	Delimiter for the spectra file. (Default ',').\n"
			<< " -db	Delimiter for the band map file. (Default ',').\n"
			<< " -do	Delimiter for the output file. (Default ',').\n"
			<< " -fr	First data row index (zero-based). (Default 0).\n"
			<< " -fc 	First data column index (zero-based). (Default 0).\n"
			<< " -dc 	Date column index (zero-based). (Default -1).\n"
			<< " -tc 	Timestamp column index (zero-based). (Default -1).\n"
			<< "    Run without arguments to use the gui.\n";
}

#include "convolver.hpp"

int main(int argc, char** argv) {

	if(argc > 1) {
		if(argc < 4) {
			usage();
			return 1;
		} else {
			std::vector<std::string> files;
			double inputScale = 1.0;
			double threshold = 0.0001;
			double shift = 0;
			std::string bandDelim = ",";
			std::string specDelim = ",";
			std::string outputDelim = ",";
			int firstRow = 0;
			int firstCol = 0;
			int dateCol = -1;
			int timeCol = -1;

			for(int i = 1; i < argc; ++i) {
				std::string arg = argv[++i];
				if(arg == "-s") {
					inputScale = atof(argv[++i]);
				} else if(arg == "-t"){
					threshold = atof(argv[++i]);
				} else if(arg == "-f") {
					shift = atof(argv[++i]);
				} else if(arg == "-ds") {
					specDelim = argv[++i];
				} else if (arg == "-db") {
					bandDelim = argv[++i];
				} else if (arg == "-do") {
					outputDelim = argv[++i];
				} else if (arg == "-fr") {
					firstRow = atoi(argv[++i]);
				} else if (arg == "-fc") {
					firstCol = atoi(argv[++i]);
				} else if (arg == "-dc") {
					dateCol = atoi(argv[++i]);
				} else if (arg == "-tc") {
					timeCol = atoi(argv[++i]);
				} else {
					files.push_back(arg);
				}
			}

			std::string bandDef = files[0];
			std::string spectra = files[1];
			std::string output = files[2];

			Convolver conv;
			DummyListener listener;
			bool running = true;
			conv.run(listener, bandDef, bandDelim, spectra, specDelim, firstRow, firstCol, dateCol, timeCol, output, outputDelim, inputScale, threshold, shift, running);
		}
	} else {
		return runWithGui(argc, argv);
	}

}
