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
#include "util.hpp"

using namespace hlrg::convolve;
using namespace hlrg::util;

namespace hlrg {
namespace convolve {

class ConvolveForm : public QDialog, public Ui::ConvolveForm, public ConvolveListener {
	Q_OBJECT
private:
	QSettings m_settings;
	std::string m_bandDefFile;
	std::string m_bandDefDelim;
	std::string m_spectraFile;
	std::vector<std::string> m_spectraList;
	std::string m_spectraDelim;
	std::string m_outputFile;
	std::string m_outputDelim;
	FileType m_outputType;
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

	std::vector<QWidget*> runWidgets;		///<! Widgets which are active in the run state.
	std::vector<QWidget*> stopWidgets;		///<! Widgets which are active in stopped state.
	std::vector<QWidget*> csvInputWidgets;	///<! Widgets which are active when the spectra file is a CSV.
	std::vector<QWidget*> csvOutputWidgets; ///<! Widgets which are active when the output file is a CSV.

	std::thread m_thread;
	bool m_running;

	/**
	 * Set the enabled state of CSV-related widgets depending on what types of file are selected.
	 */
	void updateCSVWidgets();

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
	void cboOutputTypeChanged(QString);
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
