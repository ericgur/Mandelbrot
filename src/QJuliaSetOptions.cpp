/**
 * @file QJuliaSetOptions.cpp
 * @brief Implementation of the QJuliaSetOptions dialog.
 *
 * Contains 10 curated Julia set preset constants and the logic for
 * preset selection, manual input, and auto-apply behavior.
 */

#include "pch.h"
#include "QJuliaSetOptions.h"

using namespace std;

/// @brief Curated preset Julia set complex constants.
constexpr array presets = {complex<double>(0.285, 0.01), complex<double>(-0.7269, 0.1889), complex<double>(-0.8, 0.156),       complex<double>(-0.4, 0.6),
                           complex<double>(0.0, -0.8),   complex<double>(0.39, 0.18),      complex<double>(-0.70176, -0.3842), complex<double>(-0.75, 0.11),
                           complex<double>(-0.1, 0.651), complex<double>(-0.712, 0.27015)};

QJuliaSetOptions::QJuliaSetOptions(QWidget* parent) : QDialog(parent), c(presets[0])
{
    ui.setupUi(this);

    ui.real->setValidator(new QDoubleValidator(-2, 2, 10, this));
    ui.imag->setValidator(new QDoubleValidator(-2, 2, 10, this));
    ui.real->setText(QString::number(c.real()));
    ui.imag->setText(QString::number(c.imag()));
    for (const auto& p : presets) {
        ui.presets->addItem(QString::number(p.real()) + ", " + QString::number(p.imag()));
    }
    connect(ui.presets, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &QJuliaSetOptions::onPresetChanged);
    connect(ui.buttonApply, &QPushButton::clicked, this, &QJuliaSetOptions::onApplyButtonClicked);
    connect(ui.real, &QLineEdit::textChanged, this, &QJuliaSetOptions::valueChanged);
    connect(ui.imag, &QLineEdit::textChanged, this, &QJuliaSetOptions::valueChanged);
}

void QJuliaSetOptions::onApplyButtonClicked()
{
    c = std::complex<double>(ui.real->text().toDouble(), ui.imag->text().toDouble());
    emit juliaConstantChanged(c);
}

void QJuliaSetOptions::valueChanged()
{
    c = std::complex<double>(ui.real->text().toDouble(), ui.imag->text().toDouble());
    if (ui.checkAutoApply->isChecked()) {
        emit juliaConstantChanged(c);
    }
}

/**
 * @brief Load a preset constant into the input fields.
 *
 * Uses QSignalBlocker to prevent feedback loops while updating
 * the text fields. Auto-applies the value if the checkbox is checked.
 *
 * @param index Index into the presets array.
 */
void QJuliaSetOptions::onPresetChanged(int index)
{
    if (index < 0 || index >= static_cast<int>(presets.size()))
        return;
    c = presets[index];

    QSignalBlocker blocker1(ui.real);
    QSignalBlocker blocker2(ui.imag);
    ui.real->setText(QString::number(c.real()));
    ui.imag->setText(QString::number(c.imag()));

    if (ui.checkAutoApply->isChecked()) {
        emit juliaConstantChanged(c);
    }
}
