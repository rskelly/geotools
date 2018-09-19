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

#include "ui_contrem.h"
#include "contrem.hpp"

using namespace hlrg;

class ContremForm : public QDialog, public Ui::ContremForm, public ContremListener {
	Q_OBJECT
private:
	QSettings m_settings;
	std::string m_spectraFile;
	std::string m_roiFile;
	std::string m_wlMapFile;
	std::string m_outputFile;
	std::string m_outputDriver;

	bool m_header;
	int m_wlCol;
	int m_bandCol;
	double m_minWl;
	double m_maxWl;
	bool m_popStats;
	int m_buffer;
	int m_threads;

	Contrem* m_contrem;
	Reader* m_reader;
	QDialog* m_form;
	QApplication* m_app;

	std::thread m_thread;
	bool m_running;

public:
	ContremForm(Contrem* contrem, QApplication* app);
	void setupUi(QDialog* form);
	void checkRun();

	void run();
	void cancel();

	void runState();
	void stopState();

signals:
	void started(Contrem*);
	void update(Contrem*);
	void stopped(Contrem*);
	void finished(Contrem*);

public slots:
	void txtROIChanged(QString);
	void txtWLChanged(QString);
	void txtSpectraChanged(QString);
	void txtOutputChanged(QString);
	void cboDriverChanged(QString);
	void chkHeaderChanged(bool);
	void spnWLColChanged(int);
	void spnBandColChanged(int);
	void spnMinWLChanged(double);
	void spnMaxWLChanged(double);
	void rdoPopStatsChanged(bool);
	void spnBufferChanged(int);
	void spnThreadsChanged(int);
	void btnROIClicked();
	void btnWLClicked();
	void btnSpectraClicked();
	void btnOutputClicked();
	void btnRunClicked();
	void btnCancelClicked();
	void btnHelpClicked();
	void btnCloseClicked();

	void convStarted(Contrem*);
	void convStopped(Contrem*);
	void convUpdate(Contrem*);
	void convFinished(Contrem*);

};

#endif /* SRC_UI_CONVOLVE_UI_HPP_ */
