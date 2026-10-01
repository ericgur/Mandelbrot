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

    Ui_QtMainWindow ui;                              ///< Qt Designer generated UI.
    QMandelbrotWidget* m_centralWidget;              ///< Central fractal rendering widget.
    QVector<QAction*> iterActions;                   ///< Iteration menu action list.
    QActionGroup* setTypeGroup = nullptr;            ///< Exclusive group for set type selection.
    QActionGroup* iterGroup = nullptr;               ///< Exclusive group for iteration selection.
    QJuliaSetOptions* juliaOptionsDialog = nullptr;  ///< Julia set configuration dialog.
    QSlider* _iterSlider = nullptr;                  ///< Log-scale iteration limit slider.
    QLabel* _iterLabel = nullptr;                    ///< Approximate value at the slider position.
};
