/**
 * @file QJuliaSetOptions.cpp
 * @brief Implementation of the QJuliaSetOptions dialog.
 *
 * Contains 10 curated Julia set preset constants and the logic for
 * preset selection, manual input, and auto-apply behavior.
 */

#include "pch.h"
#include <QRegularExpressionValidator>
#include "Fp128Text.h"
#include "QJuliaSetOptions.h"

namespace
{

/// @brief A preset Julia constant, as decimal text so it converts to fp128 without passing through a double.
struct JuliaPreset {
    const char* real;  ///< Real part.
    const char* imag;  ///< Imaginary part.
};

/// @brief Curated preset Julia set complex constants.
constexpr JuliaPreset presets[] = {{"0.285", "0.01"}, {"-0.7269", "0.1889"},   {"-0.8", "0.156"}, {"-0.4", "0.6"},   {"0", "-0.8"},
                                   {"0.39", "0.18"},  {"-0.70176", "-0.3842"}, {"-0.75", "0.11"}, {"-0.1", "0.651"}, {"-0.712", "0.27015"}};

}  // namespace

QJuliaSetOptions::QJuliaSetOptions(QWidget* parent) : QDialog(parent), c(QMandelbrotWidget::defaultJuliaConstant())
{
    ui.setupUi(this);

    // a regular expression rather than QDoubleValidator: it takes all 36 decimals, and keeps
    // the '.' decimal point that ParseJuliaComponent() expects whatever the system locale
    auto* validator = new QRegularExpressionValidator(JuliaComponentRegularExpression(), this);
    ui.real->setValidator(validator);
    ui.imag->setValidator(validator);
    ui.real->setText(ToDecimalText(c.real));
    ui.imag->setText(ToDecimalText(c.imag));
    for (const JuliaPreset& p : presets) {
        ui.presets->addItem(QString("%1, %2").arg(p.real, p.imag));
    }
    connect(ui.presets, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &QJuliaSetOptions::onPresetChanged);
    connect(ui.buttonApply, &QPushButton::clicked, this, &QJuliaSetOptions::onApplyButtonClicked);
    connect(ui.real, &QLineEdit::textChanged, this, &QJuliaSetOptions::valueChanged);
    connect(ui.imag, &QLineEdit::textChanged, this, &QJuliaSetOptions::valueChanged);
}

void QJuliaSetOptions::setConstant(const Complex128& constant)
{
    ui.real->setText(ToDecimalText(constant.real));
    ui.imag->setText(ToDecimalText(constant.imag));
}

bool QJuliaSetOptions::ReadFields()
{
    const std::optional<fp128_t> real = ParseJuliaComponent(ui.real->text());
    const std::optional<fp128_t> imag = ParseJuliaComponent(ui.imag->text());
    if (!real || !imag) {
        return false;
    }

    c = {*real, *imag};

    return true;
}

void QJuliaSetOptions::onApplyButtonClicked()
{
    // incomplete text, such as a lone "-", has nothing to apply
    if (ReadFields()) {
        emit juliaConstantChanged(c);
    }
}

void QJuliaSetOptions::valueChanged()
{
    if (ReadFields() && ui.checkAutoApply->isChecked()) {
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
    if (index < 0 || index >= static_cast<int>(std::size(presets))) {
        return;
    }
    c = {fp128_t(presets[index].real), fp128_t(presets[index].imag)};

    QSignalBlocker blocker1(ui.real);
    QSignalBlocker blocker2(ui.imag);
    ui.real->setText(ToDecimalText(c.real));
    ui.imag->setText(ToDecimalText(c.imag));

    if (ui.checkAutoApply->isChecked()) {
        emit juliaConstantChanged(c);
    }
}
