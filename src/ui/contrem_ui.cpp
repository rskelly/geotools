/*
 * convolve.cpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#include <QtWidgets/QDialog>
#include <QtCore/QString>
#include <QtWidgets/QFileDialog>
#include <QtCore/QDir>
#include <QtGui/QDesktopServices>
#include <QtWidgets/QMessageBox>

#include "contrem_ui.hpp"
#include "contrem.hpp"

ContremForm::ContremForm(Contrem* contrem, QApplication* app) :
	m_header(false),
	m_wlCol(0),
	m_bandCol(0),
	m_minWl(0),
	m_maxWl(0),
	m_popStats(false),
	m_buffer(256),
	m_threads(1),
	m_contrem(contrem),
	m_reader(nullptr),
	m_form(nullptr),
	m_app(app),
	m_running(false) {
}

void ContremForm::setupUi(QDialog* form) {
	Ui::ContremForm::setupUi(form);
	m_form = form;
	progressBar->setValue(0);

	QStringList drivers;
	drivers << "ENVI" << "GTiff";
	cboDriver->addItems(drivers);

	connect(txtROI, SIGNAL(textChanged(QString)), this, SLOT(txtROIChanged(QString)));
	connect(txtWavelength, SIGNAL(textChanged(QString)), this, SLOT(txtWLChanged(QString)));
	connect(txtSpectra, SIGNAL(textChanged(QString)), this, SLOT(txtSpectraChanged(QString)));
	connect(txtOutput, SIGNAL(textChanged(QString)), this, SLOT(txtOutputChanged(QString)));
	connect(chkHeader, SIGNAL(toggled(bool)), this, SLOT(chkHeaderChanged(bool)));
	connect(spnWLCol, SIGNAL(valueChanged(int)), this, SLOT(spnWLColChanged(int)));
	connect(spnBandCol, SIGNAL(valueChanged(int)), this, SLOT(spnBandColChanged(int)));
	connect(cboDriver, SIGNAL(currentTextChanged(QString)), this, SLOT(cboDriverChanged(QString)));
	connect(spnMinWL, SIGNAL(valueChanged(double)), this, SLOT(spnMinWLChanged(double)));
	connect(spnMaxWL, SIGNAL(valueChanged(double)), this, SLOT(spnMaxWLChanged(double)));
	connect(rdoPopulation, SIGNAL(toggled(bool)), this, SLOT(rdoPopStatsChanged(bool)));
	connect(btnROI, SIGNAL(clicked()), this, SLOT(btnROIClicked()));
	connect(btnWavelength, SIGNAL(clicked()), this, SLOT(btnWLClicked()));
	connect(btnSpectra, SIGNAL(clicked()), this, SLOT(btnSpectraClicked()));
	connect(btnOutput, SIGNAL(clicked()), this, SLOT(btnOutputClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(btnRunClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(btnCancelClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(btnHelpClicked()));
	connect(btnClose, SIGNAL(clicked()), this, SLOT(btnCloseClicked()));

	connect(this, SIGNAL(started(Contrem*)), this, SLOT(convStarted(Contrem*)));
	connect(this, SIGNAL(stopped(Contrem*)), this, SLOT(convStopped(Contrem*)));
	connect(this, SIGNAL(update(Contrem*)), this, SLOT(convUpdate(Contrem*)));
	connect(this, SIGNAL(finished(Contrem*)), this, SLOT(convFinished(Contrem*)));

	m_roiFile = m_settings.value("lastROI", "").toString().toStdString();
	m_wlMapFile = m_settings.value("lastWL", "").toString().toStdString();
	m_spectraFile = m_settings.value("lastSpectra", "").toString().toStdString();
	m_outputFile = m_settings.value("lastOutput", "").toString().toStdString();
	m_outputDriver = m_settings.value("lastDriver", "ENVI").toString().toStdString();
	m_wlCol = m_settings.value("lastWLCol", 0).toInt();
	m_bandCol = m_settings.value("lastBandCol", 0).toInt();
	m_header = m_settings.value("lastHeader", false).toBool();
	m_minWl = m_settings.value("lastMinWL", 0).toDouble();
	m_maxWl = m_settings.value("lastMaxWL", 0).toDouble();
	m_popStats = m_settings.value("lastPopStats", true).toBool();
	m_buffer = m_settings.value("lastBuffer", 256).toInt();
	m_threads = m_settings.value("lastThreads", 1).toInt();

	txtROI->setText(QString(m_roiFile.c_str()));
	txtWavelength->setText(QString(m_wlMapFile.c_str()));
	txtSpectra->setText(QString(m_spectraFile.c_str()));
	txtOutput->setText(QString(m_outputFile.c_str()));
	cboDriver->setCurrentText(QString(m_outputDriver.c_str()));
	spnWLCol->setValue(m_wlCol);
	spnBandCol->setValue(m_bandCol);
	chkHeader->setChecked(m_header);
	spnMinWL->setValue(m_minWl);
	spnMaxWL->setValue(m_maxWl);
	rdoPopulation->setChecked(m_popStats);
	spnBuffer->setValue(m_buffer);
	spnThreads->setValue(m_threads);
}

void ContremForm::checkRun() {
	/*
	bool a = !m_roiFile.empty() && QFile(m_roiFile.c_str()).exists();
	bool b = !m_spectraFile.empty() && QFile(m_spectraFile.c_str()).exists();
	QFileInfo dir(m_outputFile.c_str());
	bool c = !m_outputFile.empty() && dir.dir().exists();
	btnRun->setEnabled(a && b && c && d);
	*/
	btnRun->setEnabled(true);
}

void ContremForm::chkHeaderChanged(bool value) {
	m_header = value;
	m_settings.setValue("lastHeader", value);
	checkRun();
}

void ContremForm::rdoPopStatsChanged(bool value) {
	m_popStats = value;
	m_settings.setValue("lastPopStats", value);
	checkRun();
}

void ContremForm::spnThreadsChanged(int value) {
	m_threads = value;
	m_settings.setValue("lastThreads", value);
	checkRun();
}

void ContremForm::spnBufferChanged(int value) {
	m_buffer = value;
	m_settings.setValue("lastValue", value);
	checkRun();
}

void ContremForm::cboDriverChanged(QString value) {
	m_outputDriver = value.toStdString();
	m_settings.setValue("lastDriver", value);
	checkRun();
}

void ContremForm::txtROIChanged(QString filename) {
	m_roiFile = filename.toStdString();
	m_settings.setValue("lastROI", filename);
	checkRun();
}

void ContremForm::txtWLChanged(QString filename) {
	m_wlMapFile = filename.toStdString();
	m_settings.setValue("lastWL", filename);
	checkRun();
}

void ContremForm::txtSpectraChanged(QString filename) {
	m_spectraFile = filename.toStdString();
	m_settings.setValue("lastSpectra", filename);
	checkRun();
}

void ContremForm::txtOutputChanged(QString filename) {
	m_outputFile = filename.toStdString();
	m_settings.setValue("lastOutput", filename);
	checkRun();
}

void ContremForm::spnMinWLChanged(double value) {
	m_minWl = value;
	m_settings.setValue("lastMinWL", value);
	checkRun();
}

void ContremForm::spnMaxWLChanged(double value) {
	m_maxWl = value;
	m_settings.setValue("lastMaxWL", value);
	checkRun();
}

void ContremForm::spnWLColChanged(int value) {
	m_wlCol = value;
	m_settings.setValue("lastWLCol", value);
	checkRun();
}

void ContremForm::spnBandColChanged(int value) {
	m_bandCol = value;
	m_settings.setValue("lastBandCol", value);
	checkRun();
}

void ContremForm::btnROIClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "ENVI ROI File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtROI->setText(filename);
}

void ContremForm::btnWLClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "ENVI ROI File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtWavelength->setText(filename);
}

void ContremForm::btnSpectraClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "Spectra File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtSpectra->setText(filename);
}

void ContremForm::btnOutputClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getSaveFileName(this, "Output File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtOutput->setText(filename);
}

void _run(ContremListener* form, Contrem* contrem, Reader* reader) {
	contrem->run(form, reader);
}

void ContremForm::runState() {
	btnRun->setEnabled(false);
	btnCancel->setEnabled(true);
	txtROI->setEnabled(false);
	txtWavelength->setEnabled(false);
	btnROI->setEnabled(false);
	btnWavelength->setEnabled(false);
	txtSpectra->setEnabled(false);
	btnSpectra->setEnabled(false);
	txtOutput->setEnabled(false);
	btnOutput->setEnabled(false);
	spnWLCol->setEnabled(false);
	spnBandCol->setEnabled(false);
	chkHeader->setEnabled(false);
	cboDriver->setEnabled(false);
	spnMinWL->setEnabled(false);
	spnMaxWL->setEnabled(false);
	rdoPopulation->setEnabled(false);
	rdoSample->setEnabled(false);
	spnBuffer->setEnabled(false);
	spnThreads->setEnabled(false);
}

void ContremForm::stopState() {
	btnRun->setEnabled(true);
	btnCancel->setEnabled(false);
	txtROI->setEnabled(true);
	txtWavelength->setEnabled(true);
	btnROI->setEnabled(true);
	btnWavelength->setEnabled(true);
	txtSpectra->setEnabled(true);
	btnSpectra->setEnabled(true);
	txtOutput->setEnabled(true);
	btnOutput->setEnabled(true);
	spnWLCol->setEnabled(true);
	spnBandCol->setEnabled(true);
	chkHeader->setEnabled(true);
	cboDriver->setEnabled(true);
	spnMinWL->setEnabled(true);
	spnMaxWL->setEnabled(true);
	rdoPopulation->setEnabled(true);
	rdoSample->setEnabled(true);
	spnBuffer->setEnabled(true);
	spnThreads->setEnabled(true);
}

void ContremForm::run() {
	runState();
	if(!m_running) {
		m_running = true;

		m_contrem->bufferSize = m_buffer;
		m_contrem->driver = m_outputDriver;
		m_contrem->outfile = m_outputFile;
		m_contrem->sampleStats = !m_popStats;
		m_contrem->threads = m_threads;

		if(m_reader)
			delete m_reader;
		if(!m_roiFile.empty()) {
			m_reader = new ROIReader(m_roiFile);
		} else if(!m_spectraFile.empty()) {
			m_reader = new GDALReader(m_spectraFile);
		} else {
			throw std::invalid_argument("No input file (-r or -d) given.");
		}

		if(!m_wlMapFile.empty()) {
			if(m_wlCol == -1 || m_bandCol == -1 || m_wlCol == m_bandCol)
				throw std::invalid_argument("If a band file is given, the indices for "
						"wavelength and band must be >=0 and different from each other.");
			BandMapReader br(m_wlMapFile, m_wlCol, m_bandCol, m_header);
			m_reader->setBandMap(br.bandMap());
		}

		if(m_minWl > 0 && m_maxWl > 0)
			m_reader->setBandRange(m_minWl, m_maxWl);

		m_reader->setBufSize(m_buffer);

		m_thread = std::thread(_run, this, m_contrem, m_reader);
	}
	if(!m_thread.joinable()) {
		delete m_reader;
		m_reader = nullptr;
		stopState();
	}
}

void ContremForm::cancel() {
	if(m_running) {
		m_running = false;
		if(m_thread.joinable())
			m_thread.join();
		delete m_reader;
		m_reader = nullptr;
	}
	stopState();
}

void ContremForm::btnRunClicked() {
	run();
}

void ContremForm::btnCancelClicked() {
	cancel();
}

void ContremForm::btnHelpClicked() {
	QDesktopServices::openUrl(QUrl("https://github.com/rskelly/contrem/wiki/contrem", QUrl::TolerantMode));
}

void ContremForm::btnCloseClicked() {
	cancel();
	m_form->close();
	m_app->quit();
}

void ContremForm::convStarted(Contrem*) {
	progressBar->setValue(0);
}

void ContremForm::convUpdate(Contrem* conv) {
	progressBar->setValue(conv->progress() * 100);
}

void ContremForm::convStopped(Contrem*) {
	progressBar->setValue(0);
	stopState();
	checkRun();
}

void ContremForm::convFinished(Contrem*) {
	progressBar->setValue(100);
	QMessageBox::information(this, "Finished", "Processing is finished.");
	stopState();
	checkRun();
}
