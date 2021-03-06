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

#include "ui/convolve_ui.hpp"
#include "convolve.hpp"

using namespace hlrg::convolve;

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
	Convolve conv;
	ConvolveForm form(&conv, &q);
	QDialog qform;
	form.setupUi(&qform);
	qform.show();
	return q.exec();
}

class DummyListener : public ConvolveListener {
private:
	int lastP;
public:
	void started(Convolve*) {
		std::cout << "Running ";
		lastP = -1;
	}
	void update(Convolve* conv) {
		int p = (int) (conv->progress() * 100);
		if(p != lastP && p % 10 == 0)
			std::cout << " " << p << "% ";
		if(p != lastP) {
			std::cout << ".";
			std::cout.flush();
		}
		lastP = p;
	}
	void stopped(Convolve*) {
		std::cout << " Stopped.\n";
	}
	void finished(Convolve*) {
		std::cout << " Done.\n";
	}
};

void usage() {
	std::cerr << "Usage: convolve [[options] <band definition file> <output folder> <spectra file [spectra file [...]]> ]\n"
			<< " -t		The threshold -- the gaussian will extend out until it is below this value. (Default 0.0001).\n"
			<< " -f 	The wavelength shift. Used to counter mis-calibration in the spectrometer. (Default 0).\n"
			<< " -s		Scale the input. (Default 1.)\n"
			<< " -ds	Delimiter for the spectra file. (Default ','.)\n"
			<< " -db	Delimiter for the band map file. (Default ','.)\n"
			<< " -do	Delimiter for the output file. (Default ','.)\n"
			<< " -fr	First data row index (zero-based). (Default 0.)\n"
			<< " -fc 	First data column index (zero-based). (Default 0.)\n"
			<< " -dc 	Date column index (zero-based). (Default -1.)\n"
			<< " -tc 	Timestamp column index (zero-based). (Default -1.)\n"
			<< " -ot 	Output file type. 'CSV', 'ENVI' or 'GTiff'. (Default 'CSV'.) \n"
			<< " -m <m> The memory limit above which file-backed storage is used. (Default 0.)\n"
			<< " -p <p> Run using the given number of threads. The -m argument is multiplied by this number. (Default 1.)\n"
			<< " -g     Print the expected memory consumption given the input file(s).\n"
			<< "     Run without arguments to use the gui.\n";
}

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
			FileType outputType = FileType::CSV;
			int firstRow = 0;
			int firstCol = 0;
			int dateCol = -1;
			int timeCol = -1;
			bool guess = false;
			size_t memLimit = 0;
			int threads = 1;

			for(int i = 1; i < argc; ++i) {
				std::string arg = argv[i];
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
				} else if (arg == "-ot") {
					std::string type = argv[++i];
					if(type == "ENVI") {
						outputType = FileType::ENVI;
					} else if(type == "GTiff") {
						outputType = FileType::GTiff;
					} else if(type == "CSV") {
						outputType = FileType::CSV;
					} else {
						throw std::runtime_error("Invalid file type: " + type);
					}
				} else if(arg == "-m") {
					memLimit = std::stoull(argv[++i]);
				} else if(arg == "-p") {
					threads = std::stoi(argv[++i]);
					if(threads <= 0)
						threads = 1;
				} else if(arg == "-g") {
					guess = true;
				} else {
					files.push_back(arg);
				}
			}

			std::string bandDef = files[0];
			std::string output = files[1];
			std::vector<std::string> spectra(std::next(files.begin(), 2), files.end());

			Convolve conv;
			DummyListener listener;
			bool running = true;

			if(guess) {
				int m = 0;
				for(const std::string& f : spectra)
					m += conv.guess(f, specDelim);
				std::cout << m << "\n";
				return 0;
			} else {
				conv.run(listener, bandDef, bandDelim, spectra, specDelim, firstRow, firstCol, dateCol, timeCol, output, outputDelim, outputType, inputScale, threshold, shift, memLimit, threads, running);
			}
		}
	} else {
		return runWithGui(argc, argv);
	}

}
