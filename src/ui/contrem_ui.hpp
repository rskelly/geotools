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
	/*
	std::string m_spectraFile;			///<! The spectrum input. Could be a raster or spreadsheet.
	std::string m_spectraType;
	std::string m_roiFile;				///<! A ROI. Could be an ENVI ROI, a Shapefile, SQLite file or mask raster.
	std::string m_outputFile;			///<! The output file. Can be a raster or CSV.
	std::string m_outputType;			///<! The type of output file. GTiff, ENVI or CSV.

	double m_minWLCol;
	double m_maxWLCol;
	int m_wlHeaderRows;
	bool m_wlTranspose;

	double m_minWl;						///<! The minimum wavelength in the source.
	double m_maxWl;						///<! The maximum wavelength in the source.
	*/

	//int m_threads;

	Contrem m_contrem;					///<! Contrem processor object.

	QDialog* m_form;
	QApplication* m_app;

	std::thread m_thread;				///<! Processor thread.
	//bool m_running;						///<! True if currently running.

	void updateSpectraType();

	void updateOutputType();

	void updateWavelengths();

	void enableSpectraOptions(const std::string& filename);

public:

	/**
	 * Build the Contrem form.
	 *
	 * \param app The application.
	 */
	ContremForm(QApplication* app);

	/**
	 * Set up the user interface.
	 *
	 * \param form The form.
	 */
	void setupUi(QDialog* form);

	/**
	 * Check if the inputs are suitable to start processing.
	 * Throw exception otherwise.
	 */
	void checkRun();

	/**
	 * Run processing.
	 */
	void run();

	/**
	 * Cancel processing.
	 */
	void cancel();

	/**
	 * Set the UI to reflect the running state.
	 */
	void runState();

	/**
	 * Set the UI to reflect the stopped (initial) state.
	 */
	void stopState();

signals:

	/**
	 * Called when the processor has started.
	 *
	 * \param contrem The processor object.
	 */
	void started(Contrem* contrem);

	/**
	 * Called when the processor has an update.
	 *
	 * \param contrem The processor object.
	 */
	void update(Contrem* contrem);

	/**
	 * Called when the processor has stopped.
	 *
	 * \param contrem The processor object.
	 */
	void stopped(Contrem* contrem);

	/**
	 * Called when the processor has finished.
	 *
	 * \param contrem The processor object.
	 */
	void finished(Contrem* contrem);

public slots:
	void txtROIFileChanged(QString);
	void btnROIClicked();

	void spnMinWLColChanged(int);
	void spnMaxWLColChanged(int);
	void spnWLHeaderRowsChanged(int);
	void chkWLTransposeChanged(bool);

	void txtSpectraFileChanged(QString);
	void btnSpectraClicked();

	void txtOutputFileChanged(QString);
	void btnOutputClicked();
	void cboOutputTypeChanged(QString);

	void cboMinWLChanged(int);
	void cboMaxWLChanged(int);

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
