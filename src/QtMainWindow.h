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

    Ui_QtMainWindow ui;                              ///< Qt Designer generated UI.
    QMandelbrotWidget* m_centralWidget;              ///< Central fractal rendering widget.
    QVector<QAction*> iterActions;                   ///< Iteration menu action list.
    QActionGroup* setTypeGroup = nullptr;            ///< Exclusive group for set type selection.
    QActionGroup* iterGroup = nullptr;               ///< Exclusive group for iteration selection.
    QJuliaSetOptions* juliaOptionsDialog = nullptr;  ///< Julia set configuration dialog.
};
