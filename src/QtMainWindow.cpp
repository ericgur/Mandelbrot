/**
 * @file QtMainWindow.cpp
 * @brief Implementation of the QtMainWindow class.
 *
 * Sets up the main window UI, menu actions, keyboard navigation,
 * and status bar updates for the qMandelbrot application.
 */

#include "pch.h"
#include <cmath>
#include <QKeyEvent>
#include <QActionGroup>
#include "QtMainWindow.h"
#include "QMandelbrotWidget.h"

QtMainWindow::QtMainWindow(QWidget* parent) : QMainWindow(parent), m_centralWidget(new QMandelbrotWidget(this))
{
    ui.setupUi(this);
    setWindowTitle("Mandelbrot (Qt6)");
    setCentralWidget(m_centralWidget);
    juliaOptionsDialog = new QJuliaSetOptions(this);
    createActions();

    // initialize the widget
    onActionIterations();
}

QtMainWindow::~QtMainWindow() {}

/**
 * @brief Create and wire up all menu actions and action groups.
 *
 * Connects file save, view reset, palette animation, precision selection,
 * set type switching, iteration limits, Julia options, and OpenMP toggle
 * to their respective handlers.
 */
void QtMainWindow::createActions()
{
    // File -> Save Image resolutions
    connect(ui.actionSaveImage1920x1080, &QAction::triggered, this, &QtMainWindow::onActionSaveImage);
    connect(ui.actionSaveImage2560x1440, &QAction::triggered, this, &QtMainWindow::onActionSaveImage);
    connect(ui.actionSaveImage3840x2160, &QAction::triggered, this, &QtMainWindow::onActionSaveImage);

    // View actions
    connect(ui.actionResetZoom, &QAction::triggered, m_centralWidget, &QMandelbrotWidget::resetView);
    connect(ui.actionAnimatePalette, &QAction::toggled, m_centralWidget, &QMandelbrotWidget::setAnimatePalette);

    // catch the renderDone event and update status bar
    connect(m_centralWidget, &QMandelbrotWidget::renderDone, this, &QtMainWindow::onRenderDone);

    // Precision actions
    QActionGroup* precisionGroup = new QActionGroup(this);
    precisionGroup->addAction(ui.actionPrecisionAuto);
    precisionGroup->addAction(ui.actionDouble);
    precisionGroup->addAction(ui.actionFixedPoint128);
    precisionGroup->setExclusive(true);
    connect(ui.actionPrecisionAuto, &QAction::triggered, this, &QtMainWindow::onActionPrecision);
    connect(ui.actionFixedPoint128, &QAction::triggered, this, &QtMainWindow::onActionPrecision);
    connect(ui.actionDouble, &QAction::triggered, this, &QtMainWindow::onActionPrecision);

    // set type
    setTypeGroup = new QActionGroup(this);
    setTypeGroup->addAction(ui.actionTypeMandelbrot);
    setTypeGroup->addAction(ui.actionTypeJulia);
    setTypeGroup->setExclusive(true);
    connect(ui.actionTypeMandelbrot, &QAction::triggered, this, &QtMainWindow::onActionSetType);
    connect(ui.actionTypeJulia, &QAction::triggered, this, &QtMainWindow::onActionSetType);

    // Iterations actions
    iterActions = {ui.actionIterAuto, ui.actionIter128, ui.actionIter192,  ui.actionIter256,  ui.actionIter384,
                   ui.actionIter512,  ui.actionIter768, ui.actionIter1024, ui.actionIter1536, ui.actionIter2048};

    iterGroup = new QActionGroup(this);
    for (QAction* act : iterActions) {
        connect(act, &QAction::triggered, this, &QtMainWindow::onActionIterations);
        iterGroup->addAction(act);
    }
    iterGroup->setExclusive(true);

    // Julia constant options
    connect(ui.actionJuliaSetOptions, &QAction::triggered, juliaOptionsDialog, &QJuliaSetOptions::show);
    connect(juliaOptionsDialog, &QJuliaSetOptions::juliaConstantChanged, m_centralWidget, &QMandelbrotWidget::setJuliaConstant);

    // OpenMP support
    connect(ui.actionOpenMP, &QAction::toggled, m_centralWidget, &QMandelbrotWidget::setOpenMp);
    m_centralWidget->setOpenMp(ui.actionOpenMP->isChecked());

    // default precision
    ui.actionPrecisionAuto->setChecked(true);
}

/**
 * @brief Save the current fractal as a PNG image.
 *
 * Determines the target resolution from the triggering action and
 * delegates to QMandelbrotWidget::saveImage().
 */
void QtMainWindow::onActionSaveImage()
{
    QAction* act = qobject_cast<QAction*>(sender());
    if (!act)
        return;
    // forward to widget with resolution based on action text
    if (act == ui.actionSaveImage1920x1080)
        m_centralWidget->saveImage(1920, 1080);
    else if (act == ui.actionSaveImage2560x1440)
        m_centralWidget->saveImage(2560, 1440);
    else if (act == ui.actionSaveImage3840x2160)
        m_centralWidget->saveImage(3840, 2160);
}

/**
 * @brief Forward the selected precision mode to the rendering widget.
 *
 * Maps the triggered action to the corresponding Precision enum value.
 */
void QtMainWindow::onActionPrecision()
{
    QAction* act = qobject_cast<QAction*>(sender());
    if (!act)
        return;
    // forward to widget with resolution based on action text
    if (act == ui.actionPrecisionAuto)
        m_centralWidget->setPrecision(QMandelbrotWidget::Precision::Auto);
    else if (act == ui.actionDouble)
        m_centralWidget->setPrecision(QMandelbrotWidget::Precision::Double);
    else if (act == ui.actionFixedPoint128)
        m_centralWidget->setPrecision(QMandelbrotWidget::Precision::FixedPoint128);
}

/**
 * @brief Parse the selected iteration action and forward to the widget.
 *
 * If the action text parses as an integer, that value is used directly.
 * Otherwise, automatic iteration scaling is enabled (maxIter = 0).
 */
void QtMainWindow::onActionIterations()
{
    QAction* act = iterGroup->checkedAction();
    if (!act)
        return;
    // forward to widget with resolution based on action text
    bool ok = false;
    int iter = act->text().toInt(&ok);
    if (ok) {
        m_centralWidget->setMaximumIterations(iter);
    } else {
        // auto
        m_centralWidget->setMaximumIterations(0);
    }
}

/**
 * @brief Switch between Mandelbrot and Julia set rendering.
 */
void QtMainWindow::onActionSetType()
{
    QAction* act = setTypeGroup->checkedAction();
    if (!act)
        return;
    if (act == ui.actionTypeMandelbrot) {
        m_centralWidget->setSetType(QMandelbrotWidget::stMandelbrot);
    } else if (act == ui.actionTypeJulia) {
        m_centralWidget->setSetType(QMandelbrotWidget::stJulia);
    }
}

/**
 * @brief Open the Julia set options dialog pre-populated with the current constant.
 */
void QtMainWindow::onActionJuliaOptions()
{
    auto c = m_centralWidget->juliaConstant();
    juliaOptionsDialog->setConstant(c);
    juliaOptionsDialog->show();
}

/**
 * @brief Format and display render statistics in the status bar.
 *
 * For zoom levels above 2^16, the zoom is displayed in power-of-two notation
 * (e.g. "2^42"). Otherwise, a decimal value is shown.
 *
 * @param stats Frame statistics containing render time, zoom, size, and iterations.
 */
void QtMainWindow::onRenderDone(FrameStats stats)
{
    // if zoom is > 65536, show as power of 2
    QString zoomStr;
    const double threshold = std::pow(2.0, 16);
    if (stats.zoom > threshold) {
        double exp = std::log2(stats.zoom);
        int expRounded = static_cast<int>(std::round(exp));
        // if exponent is essentially integer, show as exact power, otherwise show rounded exponent with approximate value
        if (std::fabs(exp - expRounded) < 0.01) {
            zoomStr = QString("2^%1").arg(expRounded);
        } else {
            zoomStr = QString("2^%1 (~%2)").arg(expRounded).arg(stats.zoom, 0, 'e', 2);
        }
    } else {
        zoomStr = QString::number(stats.zoom, 'f', 2);
    }

    QString msg = QString("Render time: %1 ms | Zoom: x%2 | Size: %3x%4 | Iterations: %5")
                      .arg(stats.render_time_ms)
                      .arg(zoomStr)
                      .arg(stats.size.width())
                      .arg(stats.size.height())
                      .arg(stats.max_iterations);
    statusBar()->showMessage(msg);
}

/**
 * @brief Handle keyboard navigation events.
 *
 * Arrow keys pan the viewport by 5%. Plus/minus zoom in/out by 2x
 * from the center of the viewport.
 *
 * @param event The key event to process.
 */
void QtMainWindow::keyPressEvent(QKeyEvent* event)
{
    constexpr double panAmount = 0.05;
    switch (event->key()) {
    case Qt::Key_Up:
        m_centralWidget->panVertical(-panAmount);
        break;
    case Qt::Key_Down:
        m_centralWidget->panVertical(panAmount);
        break;
    case Qt::Key_Right:
        m_centralWidget->panHorizontal(panAmount);
        break;
    case Qt::Key_Left:
        m_centralWidget->panHorizontal(-panAmount);
        break;
    case Qt::Key_Plus:
        m_centralWidget->zoomIn();
        break;
    case Qt::Key_Minus:
        m_centralWidget->zoomOut();
        break;
    default:
        break;
    }

    QMainWindow::keyPressEvent(event);
}
