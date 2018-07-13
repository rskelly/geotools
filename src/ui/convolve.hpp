/*
 * convolve.hpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#ifndef SRC_UI_CONVOLVE_HPP_
#define SRC_UI_CONVOLVE_HPP_

#include <QtWidgets/QDialog>
#include <QtCore/QSettings>

#include "ui_convolve.h"
#include "convolver.hpp"

class ConvolveForm : public QDialog, public Ui::ConvolveForm, public ConvolverListener {
	Q_OBJECT
private:
	QSettings m_settings;
	std::string m_bandDefFile;
	std::string m_spectraFile;
	std::string m_outputFile;
	Convolver* m_convolver;
	QDialog* m_form;
	QApplication* m_app;

public:
	ConvolveForm(Convolver* convolver, QApplication* app);
	void setupUi(QDialog* form);
	void checkRun();

	void started(Convolver*);
	void update(Convolver*);
	void stopped(Convolver*);

public slots:
	void txtBandDefChanged(QString);
	void txtSpectraChanged(QString);
	void txtOutputChanged(QString);
	void btnBandDefClicked();
	void btnSpectraClicked();
	void btnOutputClicked();
	void btnRunClicked();
	void btnCancelClicked();
	void btnHelpClicked();
	void btnCloseClicked();

};

#endif /* SRC_UI_CONVOLVE_HPP_ */
