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
#include "ui/convolve.hpp"


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
		std::cout << "Started\n";
	}
	void update(Convolver* conv) {
		std::cout << "Progress: " << (int) (conv->progress() * 100) << "%\n";
	}
	void stopped(Convolver* conv) {
		std::cout << "Stopped.\n";
	}
	void finished(Convolver* conv) {
		std::cout << "Finished.\n";
	}
};

void usage() {
	std::cerr << "Usage: convolve [<band definition file> <spectra file> <output file>]\n"
			<< "    Run without arguments to use the gui.\n";
}

int main(int argc, char** argv) {

	if(argc > 1) {
		if(argc != 4) {
			usage();
			return 1;
		} else {
			std::string bandDef = argv[1];
			std::string spectra = argv[2];
			std::string output = argv[3];

			Convolver conv;
			DummyListener listener;
			bool running = true;
			conv.run(listener, bandDef, spectra, output, running);
		}
	} else {
		return runWithGui(argc, argv);
	}

}
