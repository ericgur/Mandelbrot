/**
 * @file QtMainWindow.h
 * @brief Main application window for the qMandelbrot fractal renderer.
 *
 * Provides the top-level window with menu bar, status bar, and keyboard
 * navigation. Routes user actions (precision, iterations, set type, save)
 * to the central QMandelbrotWidget rendering engine.
 */

#pragma once

#include <QMainWindow>
#include "ui_QtMainWindow.h"
#include "Favorites.h"
#include "QJuliaSetOptions.h"

class QAction;
class QAbstractButton;
class QLabel;
class QSlider;
class QMandelbrotWidget;

/**
 * @class QtMainWindow
 * @brief Main window managing menus, keyboard input, and status display.
 *
 * Owns the central QMandelbrotWidget and wires up all menu actions
 * (file save, precision selection, iteration limits, set type switching,
 * Julia options, OpenMP toggle). Displays render statistics in the status bar
 * after each frame and handles keyboard navigation (arrow keys for panning,
 * +/- for zooming).
 *
 * @par Iterations Slider
 * A vertical slider on the left of the fractal sets the iteration limit on a
 * logarithmic scale from 64 to 2500. The limit used for drawing comes from
 * whichever control the user touched last: moving the slider switches to that
 * value and unchecks the Iterations menu, while picking a menu item takes over
 * again. The slider follows menu presets, and while Auto is active it follows
 * the zoom-derived limit after each render. Slider positions are discrete, so
 * the label shows the value at the slider's position. That is approximate when
 * the menu or Auto set the limit; the status bar shows the exact value.
 * Dragging renders on release only, and the slider never takes keyboard focus
 * so the arrow and +/- keys keep navigating the view.
 *
 * @par Favorites
 * The Favorites menu lists saved locations; picking one moves the view there and
 * switches the set type, the iteration limit and, for Julia favorites, the Julia
 * constant to match. The Iterations menu and slider follow the restored limit as
 * if the user had picked it there. F8 (Add Current View) saves the current view at
 * once under a generated name such as "Mandelbrot 2^42 - 2026-10-01 14:32", and
 * Edit Favorites opens QFavoritesDialog to rename, change, reorder or delete
 * entries. The list is written to the application settings after every change,
 * and loaded from them at startup.
 */
class QtMainWindow : public QMainWindow
{
    Q_OBJECT
public:
    /**
     * @brief Construct the main window.
     * @param parent Optional parent widget.
     */
    explicit QtMainWindow(QWidget* parent = nullptr);

    /** @brief Destructor. */
    ~QtMainWindow();

private slots:
    /** @brief Handle image save actions at various resolutions. */
    void onActionSaveImage();

    /** @brief Handle precision mode selection (Auto/Double/FixedPoint128). */
    void onActionPrecision();

    /** @brief Handle iteration count selection from the menu. */
    void onActionIterations();

    /** @brief Handle fractal set type switching (Mandelbrot/Julia). */
    void onActionSetType();

    /** @brief Show the Julia set options dialog. */
    void onActionJuliaOptions();

    /**
     * @brief Update the status bar with render statistics.
     * @param stats Frame statistics from the most recent render pass.
     */
    void onRenderDone(FrameStats stats);

protected:
    /**
     * @brief Handle keyboard input for navigation.
     *
     * Arrow keys pan 5% of the viewport, +/- zoom in/out 2x.
     *
     * @param event Key event to process.
     */
    void keyPressEvent(QKeyEvent* event) override;

private:
    /** @brief Create and connect all menu actions and action groups. */
    void createActions();

    /**
     * @brief Build the iterations slider panel and lay it out left of the fractal widget.
     *
     * Replaces the central widget with a container holding the slider panel and
     * the QMandelbrotWidget side by side.
     */
    void CreateIterationsPanel();

    /**
     * @brief Apply an iteration limit picked with the slider.
     *
     * Unchecks the Iterations menu, since the slider now controls the limit,
     * and forwards the value to the widget.
     *
     * @param position Slider position; 0 is the minimum iteration limit.
     */
    void OnIterationsSliderChanged(int position);

    /**
     * @brief Move the slider to the position nearest an iteration limit without applying it.
     *
     * Signals are blocked, so the move does not count as a user selection.
     * Skipped while the user is dragging the handle.
     *
     * @param iterations Iteration limit to show.
     */
    void SyncIterationsSlider(int64_t iterations);

    /**
     * @brief Show the iteration limit at a slider position in the slider label.
     * @param position Slider position; 0 is the minimum iteration limit.
     */
    void UpdateIterationsLabel(int position);

    /**
     * @brief Capture the current view as a favorite with a generated description.
     * @return The view center, zoom, iteration limit, set type and Julia constant, described as e.g. "Mandelbrot 2^42 - 2026-10-01 14:32".
     */
    [[nodiscard]] Favorite CaptureFavorite() const;

    /**
     * @brief Move the view to a favorite.
     *
     * Switches the set type, and for a Julia favorite the Julia constant, before
     * moving the view, since changing either one resets the view. Also restores
     * the iteration limit, with the Iterations menu and slider following it.
     *
     * @param favorite Location to show.
     */
    void ApplyFavorite(const Favorite& favorite);

    /**
     * @brief Switch to an iteration limit, updating the Iterations menu and slider to match.
     * @param maxIterations Iteration limit, or QMandelbrotWidget::auto_iterations for Auto.
     */
    void SelectIterationLimit(int64_t maxIterations);

    /** @brief Save the current view as a new favorite (F8). */
    void AddCurrentViewToFavorites();

    /** @brief Open the favorites editor, and store the edited list if it is accepted. */
    void EditFavorites();

    /**
     * @brief Write the favorites to the settings, reporting a failure in the status bar.
     * @return True if the settings were written.
     */
    bool StoreFavorites();

    /** @brief Recreate the favorite entries below the fixed items of the Favorites menu. */
    void RebuildFavoritesMenu();

    /**
     * @brief Show a message in the status bar for a few seconds, then restore the render statistics.
     * @param notice Message to show.
     */
    void ShowStatusNotice(const QString& notice);

    Ui_QtMainWindow ui;                              ///< Qt Designer generated UI.
    QMandelbrotWidget* m_centralWidget;              ///< Central fractal rendering widget.
    QVector<QAction*> iterActions;                   ///< Iteration menu action list.
    QActionGroup* setTypeGroup = nullptr;            ///< Exclusive group for set type selection.
    QActionGroup* iterGroup = nullptr;               ///< Exclusive group for iteration selection.
    QJuliaSetOptions* juliaOptionsDialog = nullptr;  ///< Julia set configuration dialog.
    QSlider* _iterSlider = nullptr;                  ///< Log-scale iteration limit slider.
    QLabel* _iterLabel = nullptr;                    ///< Approximate value at the slider position.
    QVector<Favorite> _favorites;                    ///< Saved locations, in menu order.
    QVector<QAction*> _favoriteActions;              ///< Favorites menu entries, rebuilt whenever the list changes.
    QString _renderStatsMessage;                     ///< Last render statistics shown in the status bar.
};
