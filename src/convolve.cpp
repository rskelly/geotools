/*
 * convolve.cpp
 *
 *  Created on: Jul 10, 2018
 *      Author: rob
 */

#include "ui/convolve.hpp"

#include <QtCore/QObject>
#include <QtWidgets/QApplication>
#include <QtWidgets/QMessageBox>

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
	ConvolveForm form(&conv);
	QDialog qform;
	form.setupUi(&qform);
	qform.show();
	return q.exec();
}


int main(int argc, char** argv) {

	return runWithGui(argc, argv);

}
