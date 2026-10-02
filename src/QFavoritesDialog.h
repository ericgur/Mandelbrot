/**
 * @file QFavoritesDialog.h
 * @brief Dialog for editing the list of favorite locations.
 *
 * Lists the favorites next to a form for the selected one, and offers adding (from scratch
 * or from the current view), duplicating, deleting and reordering entries.
 */

#pragma once

#include <functional>
#include <QDialog>
#include "Favorites.h"
#include "ui_QFavoritesDialog.h"

/**
 * @class QFavoritesDialog
 * @brief Modal editor for the favorites list.
 *
 * The left side lists the favorites by description; the right side edits the selected one:
 * description, set type, center coordinates, zoom exponent, iteration limit and, for Julia
 * favorites, the Julia constant. Buttons below create an entry at the default view (New) or
 * at the current view (Add Current), duplicate, delete and reorder entries.
 *
 * @par Iteration limit
 * The limit is either fixed, in the range the main window's slider offers, or Auto. While
 * Auto is checked the limit box is disabled and shows the limit Auto picks at the entry's
 * zoom, following zoom edits, so unchecking Auto starts from the limit the view was drawn with.
 *
 * The dialog edits its own copy of the list. The caller reads it back with favorites() once
 * exec() returns Accepted, so Cancel discards every change.
 *
 * @par Validation
 * Center coordinates are entered as decimal text, the only form that carries their full
 * 128-bit precision. A validator only lets through text that CoordinateRegularExpression()
 * matches, and an edit reaches the entry as soon as the text parses. Text that is still
 * incomplete, such as a lone "-", is not stored: OK refuses to close while a field shows such
 * text, and selecting another entry discards it.
 *
 * @par Go To
 * Go To, or double-clicking an entry, emits goToFavorite() with the entry as edited so far, so
 * the main window behind the dialog shows it. Go To commits nothing: Cancel still discards the
 * list edits, though the view stays where Go To moved it.
 */
class QFavoritesDialog : public QDialog
{
    Q_OBJECT

public:
    /**
     * @brief Construct the dialog.
     * @param currentView Returns the main window's current view as a favorite; Add Current calls it.
     * @param parent Optional parent widget.
     */
    explicit QFavoritesDialog(std::function<Favorite()> currentView, QWidget* parent = nullptr);

    /**
     * @brief Load the list to edit and select its first entry.
     * @param favorites Favorites in menu order.
     */
    void setFavorites(const QVector<Favorite>& favorites);

    /**
     * @brief Get the list as edited so far.
     * @return Favorites in menu order.
     */
    [[nodiscard]] const QVector<Favorite>& favorites() const { return _favorites; }

public slots:
    /** @brief Close the dialog with Accepted, unless a field of the selected entry holds incomplete text. */
    void accept() override;

signals:
    /**
     * @brief Emitted when the user asks to see a favorite in the main window.
     * @param favorite The selected entry, including edits made so far.
     */
    void goToFavorite(const Favorite& favorite);

private:
    /**
     * @brief Get the selected entry.
     * @return The entry selected in the list, or nullptr when nothing is selected.
     */
    [[nodiscard]] Favorite* CurrentFavorite();

    /**
     * @brief Insert an entry into the list and select it.
     * @param row Position of the new entry.
     * @param favorite Entry to insert.
     */
    void InsertFavorite(int row, const Favorite& favorite);

    /**
     * @brief Select a row and show its entry in the form.
     * @param row Row to select, or -1 to clear the selection.
     */
    void SelectRow(int row);

    /** @brief Fill the form from the selected entry, or clear and disable it when nothing is selected. */
    void ShowCurrentFavorite();

    /** @brief Enable the buttons that apply to the selected entry and its position in the list. */
    void UpdateButtons();

    /** @brief Enable the Julia constant fields only for a Julia favorite. */
    void UpdateJuliaFields();

    /** @brief Show the selected entry's iteration limit, or the one Auto picks at its zoom, and enable the box for a fixed limit only. */
    void UpdateIterationFields();

    /**
     * @brief Show the zoom multiplier that a zoom exponent stands for.
     * @param log2Zoom Zoom exponent.
     */
    void UpdateZoomLabel(int log2Zoom);

    /** @brief Append an entry at the default view, then select its description for typing. */
    void AddNew();

    /** @brief Append an entry at the main window's current view. */
    void AddCurrent();

    /** @brief Insert a copy of the selected entry right after it. */
    void DuplicateCurrent();

    /** @brief Remove the selected entry and select its neighbor. */
    void DeleteCurrent();

    /**
     * @brief Move the selected entry up or down the list.
     * @param offset Rows to move by; negative moves up.
     */
    void MoveCurrent(int offset);

    /** @brief Emit goToFavorite() for the selected entry. */
    void GoToCurrent();

    /**
     * @brief Store an edited description and show it in the list.
     * @param text New description.
     */
    void OnDescriptionEdited(const QString& text);

    /** @brief Store the set type picked in the form. */
    void OnSetTypeChanged();

    /** @brief Store the center coordinate fields that hold a complete number. */
    void OnCenterEdited();

    /**
     * @brief Store an edited zoom exponent.
     * @param log2Zoom New zoom exponent.
     */
    void OnZoomChanged(int log2Zoom);

    /**
     * @brief Store an edited fixed iteration limit.
     * @param maxIterations New iteration limit.
     */
    void OnMaxIterationsChanged(int maxIterations);

    /**
     * @brief Switch the selected entry between Auto and a fixed iteration limit.
     * @param checked True for Auto; false fixes the limit at the value the box shows.
     */
    void OnAutoIterationsToggled(bool checked);

    /** @brief Store the Julia constant fields that hold a complete number. */
    void OnJuliaConstantEdited();

    Ui_QFavoritesDialog ui;                  ///< Qt Designer generated UI.
    QVector<Favorite> _favorites;            ///< The list being edited, one entry per list row.
    std::function<Favorite()> _currentView;  ///< Captures the main window's view for Add Current.
};
