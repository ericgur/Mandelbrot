/**
 * @file QFavoritesDialog.cpp
 * @brief Implementation of the QFavoritesDialog favorites editor.
 *
 * Keeps the list widget and the favorites vector row for row. Structural changes (insert,
 * delete, move) run with the list's signals blocked and refresh the form explicitly
 * afterwards, because QListWidget reports a current-row change before a removal reaches the
 * model, when its row number would index the vector wrongly.
 */

#include "pch.h"
#include <algorithm>
#include <cmath>
#include <utility>
#include <QMessageBox>
#include <QRegularExpressionValidator>
#include <QSignalBlocker>
#include "QFavoritesDialog.h"

namespace
{

constexpr int decimalFieldChars = 40;     ///< Coordinate and Julia fields fit a sign, an integer digit, the point and 37 decimals.
constexpr int maxPlainZoomExponent = 16;  ///< Above 2^16 the multiplier is shown in scientific notation, like the status bar.

}  // namespace

QFavoritesDialog::QFavoritesDialog(std::function<Favorite()> currentView, QWidget* parent) : QDialog(parent), _currentView(std::move(currentView))
{
    ui.setupUi(this);

    ui.setType->addItem("Mandelbrot", QMandelbrotWidget::stMandelbrot);
    ui.setType->addItem("Julia Set", QMandelbrotWidget::stJulia);
    ui.zoom->setRange(QMandelbrotWidget::logMinZoom, QMandelbrotWidget::logMaxZoom);
    // same range as the iterations slider in the main window
    ui.maxIterations->setRange(static_cast<int>(QMandelbrotWidget::min_fixed_iterations), static_cast<int>(QMandelbrotWidget::max_iterations));

    // regular expressions rather than QDoubleValidator: they take all 37 decimals, and keep the
    // '.' decimal point that the parse functions expect whatever the system locale
    auto* coordinateValidator = new QRegularExpressionValidator(CoordinateRegularExpression(), this);
    ui.centerX->setValidator(coordinateValidator);
    ui.centerY->setValidator(coordinateValidator);
    // same range as the Julia Set Options dialog
    auto* juliaValidator = new QRegularExpressionValidator(JuliaComponentRegularExpression(), this);
    ui.juliaReal->setValidator(juliaValidator);
    ui.juliaImag->setValidator(juliaValidator);

    const int decimalFieldWidth = ui.centerX->fontMetrics().horizontalAdvance(QString(decimalFieldChars, u'0'));
    for (QLineEdit* field : {ui.centerX, ui.centerY, ui.juliaReal, ui.juliaImag}) {
        field->setMinimumWidth(decimalFieldWidth);
    }

    connect(ui.favoritesList, &QListWidget::currentRowChanged, this, &QFavoritesDialog::ShowCurrentFavorite);
    connect(ui.favoritesList, &QListWidget::itemDoubleClicked, this, &QFavoritesDialog::GoToCurrent);

    connect(ui.buttonNew, &QPushButton::clicked, this, &QFavoritesDialog::AddNew);
    connect(ui.buttonAddCurrent, &QPushButton::clicked, this, &QFavoritesDialog::AddCurrent);
    connect(ui.buttonDuplicate, &QPushButton::clicked, this, &QFavoritesDialog::DuplicateCurrent);
    connect(ui.buttonDelete, &QPushButton::clicked, this, &QFavoritesDialog::DeleteCurrent);
    connect(ui.buttonMoveUp, &QPushButton::clicked, this, [this]() { MoveCurrent(-1); });
    connect(ui.buttonMoveDown, &QPushButton::clicked, this, [this]() { MoveCurrent(1); });
    connect(ui.buttonGoTo, &QPushButton::clicked, this, &QFavoritesDialog::GoToCurrent);
    connect(ui.buttonBox, &QDialogButtonBox::accepted, this, &QFavoritesDialog::accept);
    connect(ui.buttonBox, &QDialogButtonBox::rejected, this, &QFavoritesDialog::reject);

    // textEdited fires for typing only, so filling the form from an entry doesn't write it back
    connect(ui.description, &QLineEdit::textEdited, this, &QFavoritesDialog::OnDescriptionEdited);
    connect(ui.centerX, &QLineEdit::textEdited, this, &QFavoritesDialog::OnCenterEdited);
    connect(ui.centerY, &QLineEdit::textEdited, this, &QFavoritesDialog::OnCenterEdited);
    connect(ui.juliaReal, &QLineEdit::textEdited, this, &QFavoritesDialog::OnJuliaConstantEdited);
    connect(ui.juliaImag, &QLineEdit::textEdited, this, &QFavoritesDialog::OnJuliaConstantEdited);
    connect(ui.setType, &QComboBox::currentIndexChanged, this, &QFavoritesDialog::OnSetTypeChanged);
    connect(ui.zoom, &QSpinBox::valueChanged, this, &QFavoritesDialog::OnZoomChanged);
    connect(ui.maxIterations, &QSpinBox::valueChanged, this, &QFavoritesDialog::OnMaxIterationsChanged);
    connect(ui.autoIterations, &QCheckBox::toggled, this, &QFavoritesDialog::OnAutoIterationsToggled);

    ShowCurrentFavorite();
}

void QFavoritesDialog::setFavorites(const QVector<Favorite>& favorites)
{
    _favorites = favorites;
    {
        const QSignalBlocker blocker(ui.favoritesList);
        ui.favoritesList->clear();
        for (const Favorite& favorite : _favorites) {
            ui.favoritesList->addItem(favorite.displayName());
        }
    }
    SelectRow(_favorites.isEmpty() ? -1 : 0);
}

void QFavoritesDialog::accept()
{
    // the entry only takes complete numbers, so incomplete text in a field never reached it
    if (CurrentFavorite()) {
        const std::pair<QLineEdit*, const char*> fields[] = {{ui.centerX, "Center X"},
                                                             {ui.centerY, "Center Y"},
                                                             {ui.juliaReal, "The Julia constant's real part"},
                                                             {ui.juliaImag, "The Julia constant's imaginary part"}};
        for (const auto& [field, name] : fields) {
            if (field->isEnabled() && !field->hasAcceptableInput()) {
                QMessageBox::warning(this, windowTitle(), QString("%1 is not a complete number.").arg(name));
                field->setFocus();
                field->selectAll();
                return;
            }
        }
    }

    QDialog::accept();
}

Favorite* QFavoritesDialog::CurrentFavorite()
{
    const int row = ui.favoritesList->currentRow();
    if (row < 0 || row >= _favorites.size()) {
        return nullptr;
    }

    return &_favorites[row];
}

void QFavoritesDialog::InsertFavorite(int row, const Favorite& favorite)
{
    _favorites.insert(row, favorite);
    {
        const QSignalBlocker blocker(ui.favoritesList);
        ui.favoritesList->insertItem(row, favorite.displayName());
    }
    SelectRow(row);
}

void QFavoritesDialog::SelectRow(int row)
{
    {
        const QSignalBlocker blocker(ui.favoritesList);
        ui.favoritesList->setCurrentRow(row);
    }
    // explicit, since the row may not have changed and then no signal would have come anyway
    ShowCurrentFavorite();
}

void QFavoritesDialog::ShowCurrentFavorite()
{
    UpdateButtons();

    const Favorite* favorite = CurrentFavorite();
    ui.detailsPanel->setEnabled(favorite != nullptr);

    const QSignalBlocker setTypeBlocker(ui.setType);
    const QSignalBlocker zoomBlocker(ui.zoom);
    const QSignalBlocker autoIterationsBlocker(ui.autoIterations);
    if (!favorite) {
        ui.description->clear();
        ui.setType->setCurrentIndex(-1);
        ui.centerX->clear();
        ui.centerY->clear();
        ui.zoom->setValue(ui.zoom->minimum());
        ui.zoomMultiplier->clear();
        ui.autoIterations->setChecked(false);
        {
            const QSignalBlocker maxIterationsBlocker(ui.maxIterations);
            ui.maxIterations->setValue(ui.maxIterations->minimum());
        }
        ui.juliaReal->clear();
        ui.juliaImag->clear();
        return;
    }

    ui.description->setText(favorite->description);
    ui.setType->setCurrentIndex(ui.setType->findData(favorite->setType));
    ui.centerX->setText(ToDecimalText(favorite->centerX));
    ui.centerY->setText(ToDecimalText(favorite->centerY));
    ui.zoom->setValue(favorite->log2Zoom);
    UpdateZoomLabel(favorite->log2Zoom);
    ui.autoIterations->setChecked(favorite->maxIterations == QMandelbrotWidget::auto_iterations);
    UpdateIterationFields();
    ui.juliaReal->setText(ToDecimalText(favorite->juliaConstant.real));
    ui.juliaImag->setText(ToDecimalText(favorite->juliaConstant.imag));
    UpdateJuliaFields();
}

void QFavoritesDialog::UpdateButtons()
{
    const int row = ui.favoritesList->currentRow();
    const bool selected = row >= 0;
    ui.buttonDuplicate->setEnabled(selected);
    ui.buttonDelete->setEnabled(selected);
    ui.buttonGoTo->setEnabled(selected);
    ui.buttonMoveUp->setEnabled(selected && row > 0);
    ui.buttonMoveDown->setEnabled(selected && row < _favorites.size() - 1);
}

void QFavoritesDialog::UpdateJuliaFields()
{
    const Favorite* favorite = CurrentFavorite();
    const bool isJulia = favorite && favorite->setType == QMandelbrotWidget::stJulia;
    ui.labelJulia->setEnabled(isJulia);
    ui.labelJuliaImag->setEnabled(isJulia);
    ui.juliaReal->setEnabled(isJulia);
    ui.juliaImag->setEnabled(isJulia);
}

void QFavoritesDialog::UpdateIterationFields()
{
    const Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    // while Auto is on, the box shows the limit Auto picks at this zoom, so turning Auto off
    // starts from the limit the view was drawn with
    const bool isAuto = favorite->maxIterations == QMandelbrotWidget::auto_iterations;
    const int64_t limit = isAuto ? QMandelbrotWidget::autoIterationLimit(favorite->log2Zoom) : favorite->maxIterations;
    const QSignalBlocker blocker(ui.maxIterations);
    ui.maxIterations->setValue(static_cast<int>(limit));
    ui.maxIterations->setEnabled(!isAuto);
}

void QFavoritesDialog::UpdateZoomLabel(int log2Zoom)
{
    if (log2Zoom <= maxPlainZoomExponent) {
        ui.zoomMultiplier->setText(QString("(x%1)").arg(1ll << log2Zoom));
    } else {
        ui.zoomMultiplier->setText(QString("(x%1)").arg(std::ldexp(1.0, log2Zoom), 0, 'e', 2));
    }
}

void QFavoritesDialog::AddNew()
{
    Favorite favorite;
    favorite.description = "New favorite";
    InsertFavorite(static_cast<int>(_favorites.size()), favorite);

    ui.description->setFocus();
    ui.description->selectAll();
}

void QFavoritesDialog::AddCurrent()
{
    InsertFavorite(static_cast<int>(_favorites.size()), _currentView());
}

void QFavoritesDialog::DuplicateCurrent()
{
    const Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    Favorite copy = *favorite;
    copy.description += " (copy)";
    InsertFavorite(ui.favoritesList->currentRow() + 1, copy);
}

void QFavoritesDialog::DeleteCurrent()
{
    const int row = ui.favoritesList->currentRow();
    if (row < 0) {
        return;
    }

    _favorites.removeAt(row);
    {
        const QSignalBlocker blocker(ui.favoritesList);
        delete ui.favoritesList->takeItem(row);
    }
    // the entry that moved up into the deleted row, or the new last entry
    SelectRow(std::min(row, static_cast<int>(_favorites.size()) - 1));
}

void QFavoritesDialog::MoveCurrent(int offset)
{
    const int row = ui.favoritesList->currentRow();
    const int target = row + offset;
    if (row < 0 || target < 0 || target >= _favorites.size()) {
        return;
    }

    _favorites.move(row, target);
    {
        const QSignalBlocker blocker(ui.favoritesList);
        ui.favoritesList->insertItem(target, ui.favoritesList->takeItem(row));
    }
    SelectRow(target);
}

void QFavoritesDialog::GoToCurrent()
{
    if (const Favorite* favorite = CurrentFavorite()) {
        emit goToFavorite(*favorite);
    }
}

void QFavoritesDialog::OnDescriptionEdited(const QString& text)
{
    Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    favorite->description = text;
    ui.favoritesList->currentItem()->setText(favorite->displayName());
}

void QFavoritesDialog::OnSetTypeChanged()
{
    Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    favorite->setType = static_cast<QMandelbrotWidget::set_type_t>(ui.setType->currentData().toInt());
    UpdateJuliaFields();
}

void QFavoritesDialog::OnCenterEdited()
{
    Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    if (const std::optional<fp128_t> x = ParseCoordinate(ui.centerX->text())) {
        favorite->centerX = *x;
    }
    if (const std::optional<fp128_t> y = ParseCoordinate(ui.centerY->text())) {
        favorite->centerY = *y;
    }
}

void QFavoritesDialog::OnZoomChanged(int log2Zoom)
{
    Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    favorite->log2Zoom = log2Zoom;
    UpdateZoomLabel(log2Zoom);
    // an Auto limit follows the zoom
    UpdateIterationFields();
}

void QFavoritesDialog::OnMaxIterationsChanged(int maxIterations)
{
    Favorite* favorite = CurrentFavorite();
    // the box is disabled while Auto is on, and its value then only shows what Auto picks
    if (!favorite || favorite->maxIterations == QMandelbrotWidget::auto_iterations) {
        return;
    }

    favorite->maxIterations = maxIterations;
}

void QFavoritesDialog::OnAutoIterationsToggled(bool checked)
{
    Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    favorite->maxIterations = checked ? QMandelbrotWidget::auto_iterations : ui.maxIterations->value();
    UpdateIterationFields();
}

void QFavoritesDialog::OnJuliaConstantEdited()
{
    Favorite* favorite = CurrentFavorite();
    if (!favorite) {
        return;
    }

    if (const std::optional<fp128_t> real = ParseJuliaComponent(ui.juliaReal->text())) {
        favorite->juliaConstant.real = *real;
    }
    if (const std::optional<fp128_t> imag = ParseJuliaComponent(ui.juliaImag->text())) {
        favorite->juliaConstant.imag = *imag;
    }
}
