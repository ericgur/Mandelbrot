/**
 * @file QJuliaSetOptions.h
 * @brief Dialog for configuring Julia set parameters.
 *
 * Provides a dialog with 10 curated preset complex constants, custom
 * real/imaginary input fields, and an optional auto-apply mode for
 * real-time Julia set parameter updates.
 */

#pragma once

#include <complex>
#include <QDialog>
#include "ui_QJuliaSetOptions.h"

/**
 * @class QJuliaSetOptions
 * @brief Dialog for selecting or entering Julia set complex constant parameters.
 *
 * Offers 10 preset complex constants (classic spirals, branching structures, etc.)
 * via a combo box, plus manual real/imaginary input fields validated to the
 * range [-2, 2]. When auto-apply is enabled, changes are emitted immediately
 * via the juliaConstantChanged signal.
 */
class QJuliaSetOptions : public QDialog
{
    Q_OBJECT

public:
    /**
     * @brief Construct the Julia set options dialog.
     * @param parent Optional parent widget.
     */
    QJuliaSetOptions(QWidget* parent = nullptr);

    /** @brief Destructor. */
    ~QJuliaSetOptions() {}

    /**
     * @brief Populate the dialog fields with the given complex constant.
     * @param constant The complex constant to display in the real/imaginary fields.
     */
    void setConstant(std::complex<double> constant)
    {
        ui.real->setText(QString::number(constant.real()));
        ui.imag->setText(QString::number(constant.imag()));
    }

signals:
    /**
     * @brief Emitted when the Julia constant is changed and applied.
     * @param c The new complex constant value.
     */
    void juliaConstantChanged(std::complex<double> c);

private slots:
    /**
     * @brief Handle preset combo box selection changes.
     * @param index Index of the selected preset in the combo box.
     */
    void onPresetChanged(int index);

    /** @brief Apply the current real/imaginary field values and emit the signal. */
    void onApplyButtonClicked();

    /** @brief Called when real or imaginary input changes; auto-applies if enabled. */
    void valueChanged();

private:
    Ui_QJuliaSetOptions ui;  ///< Qt Designer generated UI.
    std::complex<double> c;  ///< Current complex constant value.
};
