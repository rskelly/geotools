/*
 * convolve.hpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#ifndef SRC_UI_CONVOLVE_UI_HPP_
#define SRC_UI_CONVOLVE_UI_HPP_

#include <thread>

#include <QtWidgets/QDialog>
#include <QtCore/QSettings>

#include "ui_convolve.h"
#include "convolve.hpp"

using namespace hlrg::convolve;

namespace hlrg {
namespace convolve {

class ConvolveForm : public QDialog, public Ui::ConvolveForm, public ConvolveListener {
	Q_OBJECT
private:
	QSettings m_settings;
	std::string m_bandDefFile;
	std::string m_bandDefDelim;
	std::string m_spectraFile;
	std::string m_spectraDelim;
	std::string m_outputFile;
	std::string m_outputDelim;
	double m_inputScale;
	double m_tolerance;
	double m_bandShift;
	int m_firstRow;
	int m_firstCol;
	int m_dateCol;
	int m_timeCol;

	Convolve* m_convolve;
	QDialog* m_form;
	QApplication* m_app;

	std::thread m_thread;
	bool m_running;

public:
	ConvolveForm(Convolve* convolve, QApplication* app);
	void setupUi(QDialog* form);
	void checkRun();

	void run();
	void cancel();

	void runState();
	void stopState();

signals:
	void started(Convolve*);
	void update(Convolve*);
	void stopped(Convolve*);
	void finished(Convolve*);

public slots:
	void txtBandDefChanged(QString);
	void cboBandDefDelimChanged(QString);
	void txtSpectraChanged(QString);
	void cboSpectraDelimChanged(QString);
	void spnFirstColChanged(int);
	void spnFirstRowChanged(int);
	void spnDateColChanged(int);
	void spnTimeColChanged(int);
	void txtOutputChanged(QString);
	void cboOutputDelimChanged(QString);
	void spnInputScaleChanged(double);
	void spnToleranceChanged(double);
	void spnBandShiftChanged(double);
	void btnBandDefClicked();
	void btnSpectraClicked();
	void btnOutputClicked();
	void btnRunClicked();
	void btnCancelClicked();
	void btnHelpClicked();
	void btnCloseClicked();

	void convStarted(Convolve*);
	void convStopped(Convolve*);
	void convUpdate(Convolve*);
	void convFinished(Convolve*);

	void handleException(const std::exception& ex);
};

} // convolve
} // hlrg

#endif /* SRC_UI_CONVOLVE_UI_HPP_ */
