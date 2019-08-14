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

#include "../../include/geo_util.hpp"
#include "contrem.hpp"

using namespace hlrg::contrem;
using namespace hlrg::reader;

namespace {

	constexpr const char* LAST_ROI = "lastROI";
	constexpr const char* LAST_SAMPLE_POINTS = "lastSamplePoints";
	constexpr const char* LAST_SAMPLE_POINTS_LAYER = "lastSamplePointsLayer";
	constexpr const char* LAST_SAMPLE_POINTS_ID_FIELD = "lastSamplePointsIDField";
	constexpr const char* LAST_SPECTRA = "lastSpectra";
	constexpr const char* LAST_INPUT_START_BAND = "lastInputStartBand";
	constexpr const char* LAST_INPUT_END_BAND = "lastInputEndBand";
	constexpr const char* LAST_INPUT_HAS_HEADER = "lastInputHasHeader";
	constexpr const char* LAST_INPUT_START_COL = "lastInputStartCol";
	constexpr const char* LAST_INPUT_END_COL = "lastInputEndCol";
	constexpr const char* LAST_OUTPUT = "lastOutput";
	constexpr const char* LAST_OUTPUT_TYPE = "lastOutputType";
	constexpr const char* LAST_WL_MIN_COL = "lastWLMinCol";
	constexpr const char* LAST_WL_MAX_COL = "lastWLMaxCol";
	constexpr const char* LAST_WL_BAND_COL = "lastWLBandCol";
	constexpr const char* LAST_WL_HEADER_ROWS = "lastWLHeaderRows";
	constexpr const char* LAST_WL_TRANSPOSE = "lastWLTranspose";
	constexpr const char* LAST_WL_ID_COL = "lastWLIDCol";
	constexpr const char* LAST_MIN_WL = "lastMinWl";
	constexpr const char* LAST_MAX_WL = "lastMaxWL";
	constexpr const char* LAST_BUFFER = "lastBuffer";
	constexpr const char* LAST_THREADS = "lastThreads";
	constexpr const char* LAST_DIR = "lastDir";
	constexpr const char* LAST_NORM_METHOD = "lastNormMethod";
	constexpr const char* LAST_PLOT_NORM = "lastPlotNorm";
	constexpr const char* LAST_PLOT_ORIG = "lastPlotOrig";
	constexpr const char* LAST_ONLY_SAMPLES = "lastOnlySamples";

	double nearestWl(double v, const std::map<int, double>& map) {
		if(map.empty())
			return 0;
		for(auto it = map.begin(); it != map.end(); ++it) {
			if(it->second > v) {
				if(it == map.begin()) {
					return it->second;
				} else {
					auto it0 = it; --it0;
					double a = it->second;
					double b = it0->second;
					if(std::abs(v - a) < std::abs(v - b)) {
						return a;
					} else {
						return b;
					}
				}
			}
		}
		return map.rbegin()->second;
	}

	bool __blockWlCombo;

	void trun(ContremListener* form, Contrem* contrem) {
		contrem->run(form);
	}

	/**
	 * Return a map containing pairs where the int is the 1-based band index,
	 * and the float is the wavelength. Attempts to load from raster metadata
	 * or table header. If these fail, will attempt to load from first column
	 * of presumably transposed table.
	 */
	std::map<int, double> loadWavelengths(const Contrem& contrem) {
		std::map<int, double> map;
		switch(getFileType(contrem.spectra)) {
		case FileType::GTiff:
		case FileType::ENVI:
			{
				GDALReader rdr(contrem.spectra);
				for(const auto& it : rdr.getBandMap())
					map[it.second] = (double) it.first / WL_SCALE;
			}
			break;
		case FileType::CSV:
			{
				CSVReader rdr(contrem.spectra, contrem.wlTranspose, contrem.wlHeaderRows, contrem.wlMinCol, contrem.wlMaxCol, contrem.wlIDCol);
				for(const auto& it : rdr.getBandMap())
					map[it.second] = (double) it.first / WL_SCALE;
			}
			break;
		case FileType::SHP:
		case FileType::ROI:
		default:
			std::cerr << "Missing file or invalid file type: " + contrem.spectra;
		}
		return map;
	}

}

ContremForm::ContremForm(QApplication* app) :
	m_form(nullptr),
	m_app(app) {
}

Contrem& ContremForm::contrem() {
	return m_contrem;
}

void ContremForm::setupUi(QDialog* form) {
	Ui::ContremForm::setupUi(form);
	m_form = form;
	progressBar->setValue(0);

	QStringList outputTypes;
	for(FileType type : OUTPUT_TYPES)
		outputTypes << fileTypeAsString(type).c_str();
	cboOutputType->addItems(outputTypes);

	QStringList normMethods;
	for(NormMethod method : NORM_METHODS)
		normMethods << normMethodAsString(method).c_str();
	cboNormMethod->addItems(normMethods);

	runWidgets = {
		txtROIFile, btnROI, txtSamplePoints, btnSamplePoints, txtSpectraFile, btnSpectra, spnWLHeaderRows, spnWLFirstCol, spnWLLastCol,
		spnWLIDCol, chkWLTranspose, txtOutputFile, cboOutputType, btnOutput, cboMinWL, cboMaxWL, cboNormMethod, chkPlotNorm, chkPlotOrig,
		btnRun, cboSamplePointsLayer, cboSamplePointsIDField, chkOnlySamples
	};

	stopWidgets = {
		btnCancel
	};

	connect(txtROIFile, SIGNAL(textChanged(QString)), this, SLOT(txtROIFileChanged(QString)));
	connect(btnROI, SIGNAL(clicked()), this, SLOT(btnROIClicked()));

	connect(txtSamplePoints, SIGNAL(textChanged(QString)), this, SLOT(txtSamplePointsChanged(QString)));
	connect(btnSamplePoints, SIGNAL(clicked()), this, SLOT(btnSamplePointsClicked()));
	connect(cboSamplePointsLayer, SIGNAL(currentTextChanged(QString)), this, SLOT(cboSamplePointsLayerChanged(QString)));
	connect(cboSamplePointsIDField, SIGNAL(currentTextChanged(QString)), this, SLOT(cboSamplePointsIDFieldChanged(QString)));

	connect(txtSpectraFile, SIGNAL(textChanged(QString)), this, SLOT(txtSpectraFileChanged(QString)));
	connect(btnSpectra, SIGNAL(clicked()), this, SLOT(btnSpectraClicked()));
	connect(spnWLHeaderRows, SIGNAL(valueChanged(int)), this, SLOT(spnWLHeaderRowsChanged(int)));
	connect(spnWLFirstCol, SIGNAL(valueChanged(int)), this, SLOT(spnMinWLColChanged(int)));
	connect(spnWLLastCol, SIGNAL(valueChanged(int)), this, SLOT(spnMaxWLColChanged(int)));
	connect(spnWLIDCol, SIGNAL(valueChanged(int)), this, SLOT(spnWLIDColChanged(int)));
	connect(chkWLTranspose, SIGNAL(toggled(bool)), this, SLOT(chkWLTransposeChanged(bool)));

	connect(txtOutputFile, SIGNAL(textChanged(QString)), this, SLOT(txtOutputFileChanged(QString)));
	connect(cboOutputType, SIGNAL(currentTextChanged(QString)), this, SLOT(cboOutputTypeChanged(QString)));
	connect(btnOutput, SIGNAL(clicked()), this, SLOT(btnOutputClicked()));
	connect(cboMinWL, SIGNAL(currentIndexChanged(int)), this, SLOT(cboMinWLChanged(int)));
	connect(cboMaxWL, SIGNAL(currentIndexChanged(int)), this, SLOT(cboMaxWLChanged(int)));

	connect(cboNormMethod, SIGNAL(currentTextChanged(QString)), this, SLOT(cboNormMethodChanged(QString)));
	connect(chkPlotOrig, SIGNAL(toggled(bool)), this, SLOT(chkPlotOrigChanged(bool)));
	connect(chkPlotNorm, SIGNAL(toggled(bool)), this, SLOT(chkPlotNormChanged(bool)));
	connect(chkOnlySamples, SIGNAL(toggled(bool)), this, SLOT(chkOnlySamplesChanged(bool)));

	connect(btnRun, SIGNAL(clicked()), this, SLOT(btnRunClicked()));
	connect(btnCancel, SIGNAL(clicked()), this, SLOT(btnCancelClicked()));
	connect(btnHelp, SIGNAL(clicked()), this, SLOT(btnHelpClicked()));
	connect(btnClose, SIGNAL(clicked()), this, SLOT(btnCloseClicked()));

	connect(this, SIGNAL(started(Contrem*)), this, SLOT(convStarted(Contrem*)));
	connect(this, SIGNAL(stopped(Contrem*)), this, SLOT(convStopped(Contrem*)));
	connect(this, SIGNAL(update(Contrem*)), this, SLOT(convUpdate(Contrem*)));
	connect(this, SIGNAL(finished(Contrem*)), this, SLOT(convFinished(Contrem*)));

	m_contrem.roi = m_settings.value(LAST_ROI, "").toString().toStdString();
	m_contrem.spectra = m_settings.value(LAST_SPECTRA, "").toString().toStdString();
	m_contrem.wlMinCol = m_settings.value(LAST_WL_MIN_COL, 1).toInt();
	m_contrem.wlMaxCol = m_settings.value(LAST_WL_MAX_COL, 1).toInt();
	m_contrem.wlHeaderRows = m_settings.value(LAST_WL_HEADER_ROWS, 1).toInt();
	m_contrem.wlTranspose = m_settings.value(LAST_WL_TRANSPOSE, 1).toBool();
	m_contrem.wlIDCol = m_settings.value(LAST_WL_ID_COL, -1).toInt();
	m_contrem.output = m_settings.value(LAST_OUTPUT, "").toString().toStdString();
	m_contrem.outputType = (FileType) m_settings.value(LAST_OUTPUT_TYPE, (int) FileType::ENVI).toInt();
	m_contrem.minWl = m_settings.value(LAST_MIN_WL, 0).toDouble();
	m_contrem.maxWl = m_settings.value(LAST_MAX_WL, 0).toDouble();
	m_contrem.threads = 1; //m_settings.value(LAST_THREADS, 1).toInt();
	m_contrem.samplePoints = m_settings.value(LAST_SAMPLE_POINTS, "").toString().toStdString();
	m_contrem.samplePointsLayer = m_settings.value(LAST_SAMPLE_POINTS_LAYER, "").toString().toStdString();
	m_contrem.samplePointsIDField = m_settings.value(LAST_SAMPLE_POINTS_ID_FIELD, "").toString().toStdString();
	m_contrem.plotNorm = m_settings.value(LAST_PLOT_NORM, false).toBool();
	m_contrem.plotOrig = m_settings.value(LAST_PLOT_ORIG, false).toBool();
	m_contrem.onlySamples = m_settings.value(LAST_ONLY_SAMPLES, false).toBool();
	m_contrem.normMethod = (NormMethod) m_settings.value(LAST_NORM_METHOD, (int) NormMethod::ConvexHull).toInt();

	txtROIFile->setText(QString(m_contrem.roi.c_str()));
	txtSamplePoints->setText(QString(m_contrem.samplePoints.c_str()));
	cboSamplePointsLayer->setCurrentText(QString(m_contrem.samplePointsLayer.c_str()));
	cboSamplePointsIDField->setCurrentText(QString(m_contrem.samplePointsIDField.c_str()));
	txtSpectraFile->setText(QString(m_contrem.spectra.c_str()));
	txtOutputFile->setText(QString(m_contrem.output.c_str()));
	cboOutputType->setCurrentText(QString(fileTypeAsString(m_contrem.outputType).c_str()));
	spnWLFirstCol->setValue(m_contrem.wlMinCol);
	spnWLLastCol->setValue(m_contrem.wlMaxCol);
	spnWLHeaderRows->setValue(m_contrem.wlHeaderRows);
	chkWLTranspose->setChecked(m_contrem.wlTranspose);
	chkPlotNorm->setChecked(m_contrem.plotNorm);
	chkPlotOrig->setChecked(m_contrem.plotOrig);
	chkOnlySamples->setChecked(m_contrem.onlySamples);
	cboNormMethod->setCurrentText(QString(normMethodAsString(m_contrem.normMethod).c_str()));
}

void ContremForm::checkRun() {
	FileType stype = getFileType(m_contrem.spectra);
	bool s = stype != FileType::CSV;
	txtROIFile->setEnabled(s);
	btnROI->setEnabled(s);
	txtSamplePoints->setEnabled(s);
	btnSamplePoints->setEnabled(s);
	spnWLFirstCol->setEnabled(!s);
	spnWLLastCol->setEnabled(!s);
	spnWLHeaderRows->setEnabled(!s);
	spnWLIDCol->setEnabled(!s);
	chkWLTranspose->setEnabled(!s);
	bool a = m_contrem.roi.empty() || isfile(m_contrem.roi);
	bool b = !m_contrem.spectra.empty() && isfile(m_contrem.spectra);
	bool d = m_contrem.samplePoints.empty() || isfile(m_contrem.samplePoints);
	bool c = !m_contrem.output.empty();
	bool e = !(stype == FileType::CSV && m_contrem.outputType != FileType::CSV);
	bool f = m_contrem.outputType == FileType::Unknown;
	btnRun->setEnabled(a && b && c && d);
	QStringList hints;
	if(!a)
		hints << "The ROI file is given, but doesn't exist.";
	if(!b)
		hints << "The spectra file is not given, or doesn't exist.";
	if(!c)
		hints << "The output file is not given.";
	if(!d)
		hints << "The sample points file is given, but doesn't exist.";
	if(!e)
		hints << "CSV input with raster output makes no sense.";
	if(!f)
		hints << "The output type is not valid.";

	lstHints->clear();
	lstHints->addItems(hints);
}

void ContremForm::spnMinWLColChanged(int col) {
	m_contrem.wlMinCol = col;
	m_settings.setValue(LAST_WL_MIN_COL, col);
	updateWavelengths();
	checkRun();
}

void ContremForm::spnMaxWLColChanged(int col) {
	m_contrem.wlMaxCol = col;
	m_settings.setValue(LAST_WL_MAX_COL, col);
	updateWavelengths();
	checkRun();
}

void ContremForm::spnWLHeaderRowsChanged(int rows) {
	m_contrem.wlHeaderRows = rows;
	m_settings.setValue(LAST_WL_HEADER_ROWS, rows);
	updateWavelengths();
	checkRun();
}

void ContremForm::spnWLIDColChanged(int col) {
	m_contrem.wlIDCol = col;
	m_settings.setValue(LAST_WL_ID_COL, col);
	checkRun();
}

void ContremForm::chkWLTransposeChanged(bool transpose) {
	m_contrem.wlTranspose = transpose;
	m_settings.setValue(LAST_WL_TRANSPOSE, transpose);
	updateWavelengths();
	checkRun();
}

void ContremForm::txtROIFileChanged(QString filename) {
	m_contrem.roi = filename.toStdString();
	m_settings.setValue(LAST_ROI, filename);
	checkRun();
}

void ContremForm::txtSamplePointsChanged(QString filename) {
	m_contrem.samplePoints = filename.toStdString();
	m_settings.setValue(LAST_SAMPLE_POINTS, filename);
	try {
		std::vector<std::string> names = PointSetReader::getLayerNames(filename.toStdString());
		QStringList lst;
		for(const std::string& n : names)
			lst << n.c_str();
		cboSamplePointsLayer->clear();
		cboSamplePointsLayer->addItems(lst);
		cboSamplePointsLayer->setCurrentText(QString(m_contrem.samplePointsLayer.c_str()));
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
	}
	checkRun();
}

void ContremForm::cboSamplePointsLayerChanged(QString layer) {
	m_contrem.samplePointsLayer = layer.toStdString();
	m_settings.setValue(LAST_SAMPLE_POINTS_LAYER, layer);
	try {
		std::vector<std::string> names = PointSetReader::getFieldNames(m_contrem.samplePoints, m_contrem.samplePointsLayer);
		QStringList lst;
		for(const std::string& n : names)
			lst << n.c_str();
		lst << "[auto]";
		cboSamplePointsIDField->clear();
		cboSamplePointsIDField->addItems(lst);
		cboSamplePointsIDField->setCurrentText(QString(m_contrem.samplePointsIDField.c_str()));
	} catch(const std::exception& ex) {
		std::cerr << ex.what() << "\n";
	}
	checkRun();
}

void ContremForm::cboSamplePointsIDFieldChanged(QString field) {
	m_contrem.samplePointsIDField = field.toStdString();
	m_settings.setValue(LAST_SAMPLE_POINTS_ID_FIELD, field);
	checkRun();
}

void ContremForm::updateWavelengths() {
	__blockWlCombo = true;
	cboMinWL->setEnabled(false);
	cboMaxWL->setEnabled(false);
	cboMinWL->clear();
	cboMaxWL->clear();

	std::map<int, double> map = loadWavelengths(m_contrem);

	int i = 0;
	for(auto it : map) {
		QString min, max;
		QList<QVariant> data;
		data << it.first << it.second;
		cboMinWL->insertItem(i, min.setNum(it.second, 'f', 3), QVariant(data));
		cboMaxWL->insertItem(i, max.setNum(it.second, 'f', 3), QVariant(data));
		++i;
	}

	QString min, max;
	min.setNum(nearestWl(m_contrem.minWl, map), 'f', 3);
	max.setNum(nearestWl(m_contrem.maxWl, map), 'f', 3);

	cboMinWL->setCurrentText(min);
	cboMaxWL->setCurrentText(max);

	if(!map.empty()) {
		cboMinWL->setEnabled(true);
		cboMaxWL->setEnabled(true);
	}
	__blockWlCombo = false;
}

void ContremForm::enableSpectraOptions(const std::string& spectra) {
	bool enable = FileType::CSV == getFileType(spectra);
	if(enable) {
		bool transpose;
		int header, minCol, maxCol, idCol;
		CSVReader::guessFileProperties(spectra, transpose, header, minCol, maxCol, idCol);
		spnWLFirstCol->setValue(minCol);
		spnWLLastCol->setValue(maxCol);
		spnWLHeaderRows->setValue(header);
		spnWLIDCol->setValue(idCol);
		chkWLTranspose->setChecked(transpose);
	}
}

void ContremForm::txtSpectraFileChanged(QString filename) {
	m_contrem.spectra = filename.toStdString();
	m_contrem.spectraType = fileTypeFromString(m_contrem.spectra);
	enableSpectraOptions(m_contrem.spectra);
	updateWavelengths();
	m_settings.setValue(LAST_SPECTRA, filename);
	checkRun();
}

void ContremForm::updateOutputType() {
	std::string fileType = fileTypeAsString(getFileType(m_contrem.output));
	cboOutputType->setCurrentText(fileType.c_str());
}

void ContremForm::txtOutputFileChanged(QString filename) {
	m_contrem.output = filename.toStdString();
	m_settings.setValue(LAST_OUTPUT, filename);
	checkRun();
}

void ContremForm::cboOutputTypeChanged(QString value) {
	m_contrem.outputType = fileTypeFromString(value.toStdString());
	m_settings.setValue(LAST_OUTPUT_TYPE, (int) m_contrem.outputType);
	checkRun();
}

void ContremForm::cboMinWLChanged(int index) {
	if(__blockWlCombo) return;
	QList<QVariant> v = cboMinWL->itemData(index).toList();
	m_contrem.minWl = v[1].toDouble();
	m_settings.setValue(LAST_MIN_WL, v[1].toDouble());
	checkRun();
}

void ContremForm::cboMaxWLChanged(int index) {
	if(__blockWlCombo) return;
	QList<QVariant> v = cboMinWL->itemData(index).toList();
	m_contrem.maxWl = v[1].toDouble();
	m_settings.setValue(LAST_MAX_WL, v[1].toDouble());
	checkRun();
}

void ContremForm::btnROIClicked() {
	QString lastDir = m_settings.value(LAST_DIR, "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "ROI/Mask File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue(LAST_DIR, dir.dir().absolutePath());
	txtROIFile->setText(filename);
}

void ContremForm::btnSamplePointsClicked() {
	QString lastDir = m_settings.value(LAST_DIR, "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "Sample Points File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue(LAST_DIR, dir.dir().absolutePath());
	txtSamplePoints->setText(filename);
}

void ContremForm::btnSpectraClicked() {
	QString lastDir = m_settings.value(LAST_DIR, "").toString();
	QString filename = QFileDialog::getOpenFileName(this, "Spectra File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue(LAST_DIR, dir.dir().absolutePath());
	txtSpectraFile->setText(filename);
}

void ContremForm::btnOutputClicked() {
	QString lastDir = m_settings.value(LAST_DIR, "").toString();
	QString filename = QFileDialog::getSaveFileName(this, "Output File", lastDir);
	QFileInfo dir(filename);
	m_settings.setValue(LAST_DIR, dir.dir().absolutePath());
	txtOutputFile->setText(filename);
}

void ContremForm::cboNormMethodChanged(QString str) {
	m_settings.setValue(LAST_NORM_METHOD, (int) normMethodFromString(str.toStdString()));
	checkRun();
}

void ContremForm::chkPlotOrigChanged(bool on) {
	m_settings.setValue(LAST_PLOT_ORIG, on);
	m_contrem.plotOrig = on;
	checkRun();
}

void ContremForm::chkOnlySamplesChanged(bool on) {
	m_settings.setValue(LAST_ONLY_SAMPLES, on);
	m_contrem.onlySamples = on;
	checkRun();
}

void ContremForm::chkPlotNormChanged(bool on) {
	m_settings.setValue(LAST_PLOT_NORM, on);
	m_contrem.plotNorm = on;
	checkRun();
}

void ContremForm::runState() {
	for(QWidget* w : runWidgets)
		w->setEnabled(false);
	for(QWidget* w : stopWidgets)
		w->setEnabled(true);
}

void ContremForm::stopState() {
	for(QWidget* w : runWidgets)
		w->setEnabled(true);
	for(QWidget* w : stopWidgets)
		w->setEnabled(false);
	checkRun();
}

void ContremForm::run() {
	runState();

	if(!m_contrem.running) {
		std::cout << "Starting...\n";
		m_contrem.running = true;
		m_thread = std::thread(trun, this, &m_contrem);
	}
	if(!m_thread.joinable()) {
		std::cout << "Failed to start.\n";
		m_contrem.running = false;
		stopState();
	}
}

void ContremForm::cancel() {
	if(m_contrem.running) {
		m_contrem.running = false;
		if(m_thread.joinable())
			m_thread.join();
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

void ContremForm::convStopped(Contrem* conv) {
	conv->running = false;
	progressBar->setValue(conv->progress() * 100);
	stopState();
	checkRun();
}

void ContremForm::convFinished(Contrem* conv) {
	conv->running = false;
	m_thread.detach();
	progressBar->setValue(conv->progress() * 100);
	QMessageBox::information(this, "Finished", "Processing is finished.");
	stopState();
	checkRun();
}
