#pragma once

#include <complex>
#include <QDialog>
#include "ui_QJuliaSetOptions.h"

class QJuliaSetOptions : public QDialog
{
    Q_OBJECT

public:
    QJuliaSetOptions(QWidget *parent = nullptr);
    ~QJuliaSetOptions() {}

    void setConstant(std::complex<double> constant) {
        ui.real->setText(QString::number(constant.real()));
        ui.imag->setText(QString::number(constant.imag()));
    }


signals:
    void juliaConstantChanged(std::complex<double>);

private slots:
    void onPresetChanged(int index);
    void onApplyButtonClicked();
    void valueChanged();

private:
    Ui_QJuliaSetOptions ui;
    std::complex<double> c;
};

