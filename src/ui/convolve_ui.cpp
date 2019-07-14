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

#include "convolve_ui.hpp"

using namespace hlrg::convolve;
using namespace hlrg::util;

namespace {

	constexpr const char* LAST_SPECTRA_DELIM = "lastSpectraDelim";
	constexpr const char* LAST_OUTPUT_DELIM = "lastOutputDelim";
	constexpr const char* LAST_OUTPUT_TYPE = "lastOutputType";
	constexpr const char* LAST_BAND_DELIM = "lastBandDelim";

	std::string parseDelim(const std::string& delim) {
		if(delim == "[tab]") {
			return "\t";
		} else if(delim == "[space]") {
			return " ";
		} else {
			return delim;
		}
	}

	std::string formatDelim(const std::string& delim) {
		if(delim == "\t") {
			return "[tab]";
		} else if(delim == " ") {
			return "[space]";
		} else {
			return delim;
		}
	}

	void runner(ConvolveForm* form, Convolve* conv,
			const std::string* bandDef, const std::string* bandDefDelim,
			const std::string* spectra, const std::string* spectraDelim,
			int firstRow, int firstCol, int dateCol, int timeCol,
			const std::string* output, const std::string* outputDelim, FileType outputType,
			double inputScale, double tolerance, double bandShift, bool* running) {
		try {
			conv->run(*form, *bandDef, *bandDefDelim, *spectra, *spectraDelim, firstRow, firstCol, dateCol, timeCol,
					*output, *outputDelim, outputType, inputScale, tolerance, bandShift, 0, *running);  // TODO: Mem Limit.
		} catch(const std::exception& ex) {
			form->handleException(ex);
		}
	}

}

ConvolveForm::ConvolveForm(Convolve* convolve, QApplication* app) :
	m_bandDefDelim(","),
	m_spectraDelim(","),
	m_outputDelim(","), m_outputType(FileType::Unknown),
	m_inputScale(1),
	m_tolerance(0.0001),
	m_bandShift(0),
	m_firstRow(0), m_firstCol(0),
	m_dateCol(-1), m_timeCol(-1),
	m_convolve(convolve),
	m_form(nullptr),
	m_app(app),
	m_running(false) {
}

void ConvolveForm::setupUi(QDialog* form) {
	Ui::ConvolveForm::setupUi(form);
	m_form = form;
	progressBar->setValue(0);

	stopWidgets = {
			txtBandDef, cboBandDefDelim, txtSpectra, cboSpectraDelim, spnFirstCol, spnFirstRow, spnDateCol, spnTimeCol,
			txtOutput, cboOutputDelim, cboOutputType, spnInputScale, spnTolerance, spnBandShift, btnBandDef, btnSpectra,
			btnOutput, btnRun
	};

	runWidgets = {
			btnCancel
	};

	csvInputWidgets = {
			cboSpectraDelim, spnFirstCol, spnFirstRow, spnDateCol, spnTimeCol
	};

	csvOutputWidgets = {
			cboOutputDelim
	};

	QStringList outputTypes;
	for(FileType type : OUTPUT_TYPES)
		outputTypes << fileTypeAsString(type).c_str();

	cboBandDefDelim->addItems({",", "[tab]", "[space]"});
	cboSpectraDelim->addItems({",", "[tab]", "[space]"});
	cboOutputDelim->addItems({",", "[tab]", "[space]"});
	cboOutputType->addItems(outputTypes);

	connect(txtBandDef, SIGNAL(textChanged(QString)), this, SLOT(txtBandDefChanged(QString)));
	connect(cboBandDefDelim, SIGNAL(currentTextChanged(QString)), this, SLOT(cboBandDefDelimChanged(QString)));
	connect(txtSpectra, SIGNAL(textChanged(QString)), this, SLOT(txtSpectraChanged(QString)));
	connect(cboSpectraDelim, SIGNAL(currentTextChanged(QString)), this, SLOT(cboSpectraDelimChanged(QString)));
	connect(spnFirstCol, SIGNAL(valueChanged(int)), this, SLOT(spnFirstColChanged(int)));
	connect(spnFirstRow, SIGNAL(valueChanged(int)), this, SLOT(spnFirstRowChanged(int)));
	connect(spnDateCol, SIGNAL(valueChanged(int)), this, SLOT(spnDateColChanged(int)));
	connect(spnTimeCol, SIGNAL(valueChanged(int)), this, SLOT(spnTimeColChanged(int)));
	connect(txtOutput, SIGNAL(textChanged(QString)), this, SLOT(txtOutputChanged(QString)));
	connect(cboOutputDelim, SIGNAL(currentTextChanged(QString)), this, SLOT(cboOutputDelimChanged(QString)));
	connect(cboOutputType, SIGNAL(currentTextChanged(QString)), this, SLOT(cboOutputTypeChanged(QString)));
	connect(spnInputScale, SIGNAL(valueChanged(double)), this, SLOT(spnInputScaleChanged(double)));
	connect(spnTolerance, SIGNAL(valueChanged(double)), this, SLOT(spnToleranceChanged(double)));
	connect(spnBandShift, SIGNAL(valueChanged(double)), this, SLOT(spnBandShiftChanged(double)));
	connect(btnBandDef, SIGNAL(clicked()), this, SLOT(btnBandDefClicked()));
	connect(btnSpectra, SIGNAL(clicked()), this, SLOT(btnSpectraClicked()));
	connect(btnOutput, SIGNAL(clicked()), this, SLOT(btnOutputClicked()));
	connect(btnRun, SIGNAL(clicked()), this, SLOT(btnRunClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(btnCancelClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(btnHelpClicked()));
	connect(btnClose, SIGNAL(clicked()), this, SLOT(btnCloseClicked()));

	connect(this, SIGNAL(started(Convolve*)), this, SLOT(convStarted(Convolve*)));
	connect(this, SIGNAL(stopped(Convolve*)), this, SLOT(convStopped(Convolve*)));
	connect(this, SIGNAL(update(Convolve*)), this, SLOT(convUpdate(Convolve*)));
	connect(this, SIGNAL(finished(Convolve*)), this, SLOT(convFinished(Convolve*)));

	txtBandDef->setText(m_settings.value("lastBandDef", "").toString());
	cboBandDefDelim->setCurrentText(formatDelim(m_settings.value(LAST_BAND_DELIM, ",").toString().toStdString()).c_str());
	txtSpectra->setText(m_settings.value("lastSpectra", "").toString());
	cboSpectraDelim->setCurrentText(formatDelim(m_settings.value(LAST_SPECTRA_DELIM, ",").toString().toStdString()).c_str());
	spnFirstCol->setValue(m_settings.value("lastFirstCol", 0).toInt());
	spnFirstRow->setValue(m_settings.value("lastFirstRow", 0).toInt());
	spnDateCol->setValue(m_settings.value("lastDateCol", -1).toInt());
	spnTimeCol->setValue(m_settings.value("lastTimeCol", -1).toInt());
	txtOutput->setText(m_settings.value("lastOutput", "").toString());
	cboOutputDelim->setCurrentText(formatDelim(m_settings.value(LAST_OUTPUT_DELIM, ",").toString().toStdString()).c_str());
	cboOutputType->setCurrentText(fileTypeAsString((FileType) m_settings.value(LAST_OUTPUT_TYPE, (int) FileType::Unknown).toInt()).c_str());
	spnInputScale->setValue(m_settings.value("lastInputScale", 1.0).toDouble());
	spnTolerance->setValue(m_settings.value("lastTolerance", 0.0001).toDouble());
	spnBandShift->setValue(m_settings.value("lastBandShift", 0).toDouble());

}

void ConvolveForm::updateCSVWidgets() {
	bool e = FileType::CSV == getFileType(m_spectraFile);
	for(QWidget* w : csvInputWidgets)
		w->setEnabled(e);
	bool f = FileType::CSV == getFileType(m_outputFile);
	for(QWidget* w : csvOutputWidgets)
		w->setEnabled(f);
}

void ConvolveForm::checkRun() {
	updateCSVWidgets();
	bool a = !m_bandDefFile.empty() && QFile(m_bandDefFile.c_str()).exists();
	bool b = !m_spectraFile.empty() && QFile(m_spectraFile.c_str()).exists();
	QFileInfo dir(m_outputFile.c_str());
	bool c = !m_outputFile.empty() && dir.dir().exists();
	bool d = m_tolerance > 0 && m_tolerance < 1.0;
	btnRun->setEnabled(a && b && c && d);
}

void ConvolveForm::txtBandDefChanged(QString filename) {
	m_bandDefFile = filename.toStdString();
	m_settings.setValue("lastBandDef", filename);
	checkRun();
}

void ConvolveForm::cboBandDefDelimChanged(QString delim) {
	m_bandDefDelim = parseDelim(delim.toStdString());
	m_settings.setValue(LAST_BAND_DELIM, delim);
	checkRun();
}

void ConvolveForm::txtSpectraChanged(QString filename) {
	m_spectraFile = filename.toStdString();
	m_settings.setValue("lastSpectra", filename);
	checkRun();
}

void ConvolveForm::cboSpectraDelimChanged(QString delim) {
	m_spectraDelim = parseDelim(delim.toStdString());
	m_settings.setValue(LAST_SPECTRA_DELIM, delim);
	checkRun();
}

void ConvolveForm::spnFirstColChanged(int v) {
	m_firstCol = v;
	m_settings.setValue("lastFirstCol", v);
	checkRun();
}

void ConvolveForm::spnFirstRowChanged(int v) {
	m_firstRow = v;
	m_settings.setValue("lastFirstRow", v);
	checkRun();
}

void ConvolveForm::spnDateColChanged(int v) {
	m_dateCol = v;
	m_settings.setValue("lastDateCol", v);
	checkRun();
}

void ConvolveForm::spnTimeColChanged(int v) {
	m_timeCol = v;
	m_settings.setValue("lastTimeCol", v);
	checkRun();
}

void ConvolveForm::txtOutputChanged(QString filename) {
	m_outputFile = parseDelim(filename.toStdString());
	m_settings.setValue("lastOutput", filename);
	checkRun();
}

void ConvolveForm::cboOutputDelimChanged(QString delim) {
	m_outputDelim = delim.toStdString();
	m_settings.setValue(LAST_OUTPUT_DELIM, delim);
	checkRun();
}

void ConvolveForm::cboOutputTypeChanged(QString type) {
	m_outputType = fileTypeFromString(type.toStdString());
	m_settings.setValue(LAST_OUTPUT_TYPE, (int) m_outputType);
	checkRun();
}

void ConvolveForm::spnInputScaleChanged(double value) {
	m_inputScale = value;
	m_settings.setValue("lastInputScale", value);
	checkRun();
}

void ConvolveForm::spnToleranceChanged(double value) {
	m_tolerance = value;
	m_settings.setValue("lastTolerance", value);
	checkRun();
}

void ConvolveForm::spnBandShiftChanged(double value) {
	m_bandShift = value;
	m_settings.setValue("lastBandShift", value);
	checkRun();
}

void ConvolveForm::btnBandDefClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "Band Definition File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtBandDef->setText(filename);
}

void ConvolveForm::btnSpectraClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "Spectra File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtSpectra->setText(filename);
}

void ConvolveForm::btnOutputClicked() {
	QString lastDir = m_settings.value("lastDir", "").toString();
	QString filename = QFileDialog::getSaveFileName(this, "Output File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue("lastDir", dir.dir().absolutePath());
	txtOutput->setText(filename);
}

void ConvolveForm::runState() {
	for(QWidget* w : stopWidgets)
		w->setEnabled(false);
	for(QWidget* w : runWidgets)
		w->setEnabled(true);
}

void ConvolveForm::stopState() {
	for(QWidget* w : stopWidgets)
		w->setEnabled(true);
	for(QWidget* w : runWidgets)
		w->setEnabled(false);
	updateCSVWidgets();
}

void ConvolveForm::run() {
	runState();
	if(!m_running) {
		m_running = true;
		m_thread = std::thread(runner, this, m_convolve,
				&m_bandDefFile, &m_bandDefDelim,
				&m_spectraFile, &m_spectraDelim,
				m_firstRow, m_firstCol, m_dateCol, m_timeCol,
				&m_outputFile, &m_outputDelim, m_outputType,
				m_inputScale, m_tolerance, m_bandShift, &m_running);
	}
	if(!m_thread.joinable())
		stopState();
}

bool _failure = false;

void ConvolveForm::cancel() {
	if(m_running) {
		m_running = false;
		if(_failure) {
			m_thread.detach();
		} else if(m_thread.joinable()) {
			m_thread.join();
		}
	}
	_failure = false;
	stopState();
}

void ConvolveForm::btnRunClicked() {
	run();
}

void ConvolveForm::btnCancelClicked() {
	cancel();
}

void ConvolveForm::btnHelpClicked() {
	QDesktopServices::openUrl(QUrl("https://github.com/rskelly/contrem/wiki/convolve", QUrl::TolerantMode));
}

void ConvolveForm::btnCloseClicked() {
	cancel();
	m_form->close();
	m_app->quit();
}

void ConvolveForm::convStarted(Convolve*) {
	progressBar->setValue(0);
}

void ConvolveForm::convUpdate(Convolve* conv) {
	progressBar->setValue(conv->progress() * 100);
}

void ConvolveForm::convStopped(Convolve*) {
	progressBar->setValue(0);
	stopState();
	checkRun();
}

void ConvolveForm::convFinished(Convolve*) {
	progressBar->setValue(100);
	QMessageBox::information(this, "Finished", "Processing is finished.");
	stopState();
	checkRun();
}

void ConvolveForm::handleException(const std::exception& ex) {
	progressBar->setValue(0);
	QMessageBox::critical(this, "Error", ex.what());
	_failure = true;
	cancel();
}
