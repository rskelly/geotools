/*
 * convolve.hpp
 *
 *  Created on: Jul 12, 2018
 *      Author: rob
 */

#ifndef SRC_UI_REFLECTANCE_UI_HPP_
#define SRC_UI_REFLECTANCE_UI_HPP_

#include <thread>

#include <QtWidgets/QDialog>
#include <QtCore/QSettings>

#include "../../include/reflectance.hpp"
#include "ui_reflectance.h"


namespace hlrg {

class ReflectanceForm : public QDialog, public Ui::ReflectanceForm {
	Q_OBJECT
private:
	QSettings m_settings;

	std::string m_imuGps;
	double m_imuUTCOffset;
	std::string m_rawRad;
	std::string m_frameIdx;
	std::string m_irradConv;
	double m_irradUTCOffset;
	std::string m_reflOut;

	Reflectance* m_ts;
	QDialog* m_form;
	QApplication* m_app;

	std::thread m_thread;
	bool m_running;

public:
	ReflectanceForm(Reflectance* convolver, QApplication* app);
	void setupUi(QDialog* form);
	void checkRun();

	void run();
	void cancel();

	void runState();
	void stopState();

signals:
	//void started(Convolver*);
	//void update(Convolver*);
	//void stopped(Convolver*);
	//void finished(Convolver*);

public slots:
	void txtIMUGPSChanged(QString);
	void btnIMUGPSClicked();
	void spnIMUUTCOffsetChanged(double);
	void txtFrameIndexChanged(QString);
	void btnFrameIndexClicked();
	void txtRadRastChanged(QString);
	void btnRadRastClicked();
	void txtConvIrradChanged(QString);
	void btnConvIrradClicked();
	void spnIrradUTCOffsetChanged(double);
	void txtReflOutputChanged(QString);
	void btnReflOutputClicked();
	void btnRunClicked();
	void btnCancelClicked();
	void btnHelpClicked();
	void btnCloseClicked();

	//void convStarted(Convolver*);
	//void convStopped(Convolver*);
	//void convUpdate(Convolver*);
	//void convFinished(Convolver*);

};

} // hlrg

#endif /* SRC_UI_REFLECTANCE_UI_HPP_ */