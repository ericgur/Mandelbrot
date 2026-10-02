/**
 * @file QJuliaSetOptions.h
 * @brief Dialog for configuring Julia set parameters.
 *
 * Provides a dialog with 10 curated preset complex constants, custom
 * real/imaginary input fields, and an optional auto-apply mode for
 * real-time Julia set parameter updates.
 */

#pragma once

#include <QDialog>
#include "QMandelbrotWidget.h"
#include "ui_QJuliaSetOptions.h"

/**
 * @class QJuliaSetOptions
 * @brief Dialog for selecting or entering Julia set complex constant parameters.
 *
 * Offers 10 preset complex constants (classic spirals, branching structures, etc.)
 * via a combo box, plus manual real/imaginary input fields validated to the
 * range [-2, 2]. The fields take decimal text at full fp128 precision, about 37
 * decimal places. When auto-apply is enabled, changes are emitted immediately
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
    void setConstant(const Complex128& constant);

signals:
    /**
     * @brief Emitted when the Julia constant is changed and applied.
     * @param c The new complex constant value.
     */
    void juliaConstantChanged(const Complex128& c);

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
    /**
     * @brief Read the constant from the real/imaginary fields.
     * @return True if both fields hold a complete number, which is then stored in c.
     */
    bool ReadFields();

    Ui_QJuliaSetOptions ui;  ///< Qt Designer generated UI.
    Complex128 c;            ///< Current complex constant value.
};
