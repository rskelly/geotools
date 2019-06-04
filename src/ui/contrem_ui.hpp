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
	std::string m_spectraFile;			///<! The spectrum input. Could be a raster or spreadsheet.
	std::string m_spectraType;
	std::string m_roiFile;				///<! A ROI. Could be an ENVI ROI, a Shapefile, SQLite file or mask raster.
	std::string m_roiType;				///<! A ROI. Could be an ENVI ROI, a Shapefile, SQLite file or mask raster.
	std::string m_outputFile;			///<! The output file. Can be a raster or CSV.
	std::string m_outputType;			///<! The type of output file. GTiff, ENVI or CSV.

	bool m_inputHasHeader;				///<! True if the input spreadsheet has a header that should be skipped.

	double m_minWl;						///<! The minimum wavelength in the source.
	double m_maxWl;						///<! The maximum wavelength in the source.

	int m_buffer;
	int m_threads;

	Contrem* m_contrem;					///<! Poniter to Contrem processor object.
	Reader* m_reader;					///<! Pointer to the spectral data reader.

	QDialog* m_form;
	QApplication* m_app;

	std::thread m_thread;				///<! Processor thread.
	bool m_running;						///<! True if currently running.

	void updateSpectraType();

	void updateROIType();

	void updateOutputType();

	void updateWavelengths();

public:

	/**
	 * Build the Contrem form.
	 *
	 * \param contrem The processor object.
	 * \param app The application.
	 */
	ContremForm(Contrem* contrem, QApplication* app);

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
	void cboROITypeChanged(QString);
	void btnROIClicked();

	void txtSpectraFileChanged(QString);
	void cboSpectraTypeChanged(QString);
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
