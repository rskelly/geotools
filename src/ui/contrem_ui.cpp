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

#include <gdal_priv.h>

#include "contrem_ui.hpp"
#include "contrem.hpp"

namespace {

	constexpr const char* LAST_ROI = "lastROI";
	constexpr const char* LAST_SPECTRA = "lastSpectra";
	constexpr const char* LAST_INPUT_START_BAND = "lastInputStartBand";
	constexpr const char* LAST_INPUT_END_BAND = "lastInputEndBand";
	constexpr const char* LAST_INPUT_HAS_HEADER = "lastInputHasHeader";
	constexpr const char* LAST_INPUT_START_COL = "lastInputStartCol";
	constexpr const char* LAST_INPUT_END_COL = "lastInputEndCol";
	constexpr const char* LAST_OUTPUT = "lastOutput";
	constexpr const char* LAST_OUTPUT_TYPE = "lastOutputType";
	constexpr const char* LAST_WL_COL = "lastWLCol";
	constexpr const char* LAST_WL_BAND_COL = "lastWLBandCol";
	constexpr const char* LAST_WL_HAS_HEADER = "lastWLHasHeader";
	constexpr const char* LAST_MIN_WL = "lastMinWl";
	constexpr const char* LAST_MAX_WL = "lastMaxWL";
	constexpr const char* LAST_BUFFER = "lastBuffer";
	constexpr const char* LAST_THREADS = "lastThreads";
	constexpr const char* LAST_DIR = "lastDir";

	enum FileType {
		GTiff,
		ENVI,
		ROI,
		SHP,
		CSV,
		SQLITE,
		Unknown
	};

	constexpr std::array<FileType, 3> INPUT_TYPES = {GTiff, ENVI, CSV};				///<! Allowed input types for spectral data.
	constexpr std::array<FileType, 3> OUTPUT_TYPES = {GTiff, ENVI, CSV};			///<! Allowed output types for results.
	constexpr std::array<FileType, 5> ROI_TYPES = {GTiff, ENVI, SHP, SQLITE, ROI};	///<! Allowed mask/ROI file types.

	FileType getFileType(const std::string& filename) {
		std::string ext;
		{
			size_t p = filename.find('.');
			if(p < std::string::npos) {
				std::string ext0 = filename.substr(p);
				std::transform(ext0.begin(), ext0.end(), ext.begin(), ::tolower);
			}
		}
		if(ext == ".csv") {
			return CSV;
		} else if(ext == ".roi") {
			return ROI;
		} else {
			GDALAllRegister();
			GDALDataset* ds = static_cast<GDALDataset*>(GDALOpenEx(filename.c_str(), GDAL_OF_READONLY, 0, 0, 0));
			if(ds) {
				std::string drv(ds->GetDriverName());
				FileType type = Unknown;
				if(drv == "GTiff") {
					type = GTiff;
				} else if(drv == "ENVI") {
					type = ENVI;
				} else if(drv == "ESRI Shapefile") {
					type = SHP;
				} else if(drv == "SQLite") {
					type = SQLITE;
				}
				GDALClose(ds);
				return type;
			}
		}
		return Unknown;
	}

	void trun(ContremListener* form, Contrem* contrem, Reader* reader) {
		contrem->run(form, reader);
	}

	/**
	 * Return a map containing pairs where the int is the 1-based band index,
	 * and the float is the wavelength. Attempts to load from raster metadata
	 * or table header. If these fail, will attempt to load from first column
	 * of presumably transposed table.
	 */
	std::map<int, double> loadWavelengths(const std::string& filename) {
		std::map<int, double> map;
		switch(getFileType(filename)) {
		case GTiff:
		case ENVI:
			{
				GDALReader rdr(filename);
				for(const auto& it : rdr.getBandMap())
					map[it.second] = (double) it.first / WL_SCALE;
			}
			break;
		case CSV:
			{
				CSVReader rdr(filename);
				if(true)
					rdr.transpose();
				std::vector<double> row;
				std::map<int, int> map;
				for(size_t i = 0; i < row.size(); ++i) {
					if(!std::isnan(row[i]))
						map[i] = row[i];
				}
			}
			break;
		case SHP:
		case ROI:
		default:
			throw std::runtime_error("Invalid file type: " + filename);
		}
		return map;
	}

	std::string fileTypeAsString(FileType type) {
		switch(type) {
		case GTiff: return "GTiff";
		case ENVI: return "ENVI";
		case ROI: return "ENVI ROI";
		case SHP: return "Shapefile";
		case CSV: return "CSV";
		default: return "";
		}
	}
}


ContremForm::ContremForm(Contrem* contrem, QApplication* app) :
	m_inputHasHeader(false),
	m_minWl(0),
	m_maxWl(0),
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

	QStringList outputTypes;
	outputTypes << "";
	for(FileType type : OUTPUT_TYPES)
		outputTypes << fileTypeAsString(type).c_str();
	cboOutputType->addItems(outputTypes);

	QStringList spectraTypes;
	outputTypes << "";
	for(FileType type : INPUT_TYPES)
		spectraTypes << fileTypeAsString(type).c_str();
	cboSpectraType->addItems(spectraTypes);

	QStringList roiTypes;
	roiTypes << "";
	for(FileType type : ROI_TYPES)
		roiTypes << fileTypeAsString(type).c_str();
	cboROIType->addItems(roiTypes);

	connect(txtROIFile, SIGNAL(textChanged(QString)), this, SLOT(txtROIFileChanged(QString)));
	connect(cboROIType, SIGNAL(currentTextChanged(QString)), this, SLOT(cboROITypeChanged(QString)));
	connect(txtSpectraFile, SIGNAL(textChanged(QString)), this, SLOT(txtSpectraFileChanged(QString)));
	connect(cboSpectraType, SIGNAL(currentTextChanged(QString)), this, SLOT(cboSpectraTypeChanged(QString)));
	connect(txtOutputFile, SIGNAL(textChanged(QString)), this, SLOT(txtOutputFileChanged(QString)));
	connect(cboOutputType, SIGNAL(currentTextChanged(QString)), this, SLOT(cboOutputTypeChanged(QString)));
	connect(cboMinWL, SIGNAL(currentIndexChanged(int)), this, SLOT(cboMinWLChanged(int)));
	connect(cboMaxWL, SIGNAL(currentIndexChanged(int)), this, SLOT(cboMaxWLChanged(int)));
	connect(btnROI, SIGNAL(clicked()), this, SLOT(btnROIClicked()));
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

	m_roiFile = m_settings.value(LAST_ROI, "").toString().toStdString();
	m_spectraFile = m_settings.value(LAST_SPECTRA, "").toString().toStdString();
	m_inputHasHeader = m_settings.value(LAST_INPUT_HAS_HEADER, false).toBool();
	m_outputFile = m_settings.value(LAST_OUTPUT, "").toString().toStdString();
	m_outputType = m_settings.value(LAST_OUTPUT_TYPE, "ENVI").toString().toStdString();
	m_minWl = m_settings.value(LAST_MIN_WL, 0).toDouble();
	m_maxWl = m_settings.value(LAST_MAX_WL, 0).toDouble();
	m_buffer = m_settings.value(LAST_BUFFER, 256).toInt();
	m_threads = m_settings.value(LAST_THREADS, 1).toInt();

	txtROIFile->setText(QString(m_roiFile.c_str()));
	txtSpectraFile->setText(QString(m_spectraFile.c_str()));
	txtOutputFile->setText(QString(m_outputFile.c_str()));
	cboOutputType->setCurrentText(QString(m_outputType.c_str()));
	chkInputHasHeader->setChecked(m_inputHasHeader);
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

void ContremForm::updateROIType() {
	std::string fileType = fileTypeAsString(getFileType(m_roiFile));
	cboROIType->setCurrentText(fileType.c_str());
}

void ContremForm::txtROIFileChanged(QString filename) {
	m_roiFile = filename.toStdString();
	m_settings.setValue(LAST_ROI, filename);
	updateROIType();
	checkRun();
}

void ContremForm::cboROITypeChanged(QString value) {
	m_roiType = value.toStdString();
	checkRun();
}

void ContremForm::updateSpectraType() {
	std::string fileType = fileTypeAsString(getFileType(m_spectraFile));
	cboSpectraType->setCurrentText(fileType.c_str());
}

double nearestWl(double v, const std::map<int, double>& map) {
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

void ContremForm::updateWavelengths() {
	__blockWlCombo = true;
	std::map<int, double> map = loadWavelengths(m_spectraFile);
	int i = 0;
	cboMinWL->clear();
	cboMaxWL->clear();
	for(auto it : map) {
		QString min, max;
		QList<QVariant> data;
		data << it.first << it.second;
		cboMinWL->insertItem(i, min.setNum(it.second, 'f', 3), QVariant(data));
		cboMaxWL->insertItem(i, max.setNum(it.second, 'f', 3), QVariant(data));
		++i;
	}
	QString min, max;
	min.setNum(nearestWl(m_minWl, map), 'f', 3);
	max.setNum(nearestWl(m_maxWl, map), 'f', 3);
	cboMinWL->setCurrentText(min);
	cboMaxWL->setCurrentText(max);
	__blockWlCombo = false;
}

void ContremForm::txtSpectraFileChanged(QString filename) {
	m_spectraFile = filename.toStdString();
	updateSpectraType();
	updateWavelengths();
	m_settings.setValue(LAST_SPECTRA, filename);
	checkRun();
}

void ContremForm::cboSpectraTypeChanged(QString value) {
	m_spectraType = value.toStdString();
	checkRun();
}

void ContremForm::updateOutputType() {
	std::string fileType = fileTypeAsString(getFileType(m_outputFile));
	cboOutputType->setCurrentText(fileType.c_str());
}

void ContremForm::txtOutputFileChanged(QString filename) {
	m_outputFile = filename.toStdString();
	m_settings.setValue(LAST_OUTPUT, filename);
	checkRun();
}

void ContremForm::cboOutputTypeChanged(QString value) {
	m_outputType = value.toStdString();
	checkRun();
}

void ContremForm::cboMinWLChanged(int index) {
	if(__blockWlCombo) return;
	QList<QVariant> v = cboMinWL->itemData(index).toList();
	m_minWl = v[1].toDouble();
	m_settings.setValue(LAST_MIN_WL, v[1].toDouble());
	checkRun();
}

void ContremForm::cboMaxWLChanged(int index) {
	if(__blockWlCombo) return;
	QList<QVariant> v = cboMinWL->itemData(index).toList();
	m_maxWl = v[1].toDouble();
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

void ContremForm::runState() {
	btnRun->setEnabled(false);
	btnCancel->setEnabled(true);
	txtROIFile->setEnabled(false);
	btnROI->setEnabled(false);
	txtSpectraFile->setEnabled(false);
	btnSpectra->setEnabled(false);
	txtOutputFile->setEnabled(false);
	btnOutput->setEnabled(false);
	cboOutputType->setEnabled(false);
}

void ContremForm::stopState() {
	btnRun->setEnabled(true);
	btnCancel->setEnabled(false);
	txtROIFile->setEnabled(true);
	btnROI->setEnabled(true);
	txtSpectraFile->setEnabled(true);
	btnSpectra->setEnabled(true);
	txtOutputFile->setEnabled(true);
	btnOutput->setEnabled(true);
	cboOutputType->setEnabled(true);
}

void ContremForm::run() {
	runState();
	if(!m_running) {
		m_running = true;

		m_contrem->bufferSize = m_buffer;
		m_contrem->driver = m_outputType;
		m_contrem->outfile = m_outputFile;
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

		if(m_minWl > 0 && m_maxWl > 0)
			m_reader->setBandRange(m_minWl, m_maxWl);

		m_reader->setBufSize(m_buffer);

		m_thread = std::thread(trun, this, m_contrem, m_reader);
	}
	if(!m_thread.joinable()) {
		m_running = false;
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
