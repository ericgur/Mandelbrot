/**
 * @file QtMainWindow.cpp
 * @brief Implementation of the QtMainWindow class.
 *
 * Sets up the main window UI, menu actions, keyboard navigation,
 * and status bar updates for the qMandelbrot application.
 */

#include "pch.h"
#include <algorithm>
#include <cmath>
#include <QKeyEvent>
#include <QActionGroup>
#include <QBoxLayout>
#include <QLabel>
#include <QSignalBlocker>
#include <QSlider>
#include "QtMainWindow.h"
#include "QMandelbrotWidget.h"

namespace
{

constexpr int64_t sliderMinIterations = 64;                                 ///< Iteration limit at the bottom of the slider.
constexpr int64_t sliderMaxIterations = QMandelbrotWidget::max_iterations;  ///< Iteration limit at the top of the slider.
constexpr int sliderSteps = 1000;                                           ///< Slider positions above 0; each step is ~0.37% on the log scale.
constexpr int sliderSingleStep = 5;                                         ///< Wheel step in positions (~1.9%); Qt scrolls wheelScrollLines() steps per notch.
constexpr int sliderPageStep = 50;                                          ///< Page step for groove clicks and Ctrl/Shift+wheel (~20%).

/**
 * @brief Map a slider position to an iteration limit on a logarithmic scale.
 * @param position Slider position in [0, sliderSteps].
 * @return Iteration limit in [sliderMinIterations, sliderMaxIterations].
 */
[[nodiscard]] int64_t SliderPosToIterations(int position)
{
    const double ratio = static_cast<double>(sliderMaxIterations) / sliderMinIterations;
    const double iterations = sliderMinIterations * std::pow(ratio, static_cast<double>(position) / sliderSteps);
    return std::clamp<int64_t>(std::llround(iterations), sliderMinIterations, sliderMaxIterations);
}

/**
 * @brief Map an iteration limit to the nearest slider position; the inverse of SliderPosToIterations().
 * @param iterations Iteration limit; values outside the slider range are clamped.
 * @return Slider position in [0, sliderSteps].
 */
[[nodiscard]] int IterationsToSliderPos(int64_t iterations)
{
    iterations = std::clamp(iterations, sliderMinIterations, sliderMaxIterations);
    const double ratio = static_cast<double>(sliderMaxIterations) / sliderMinIterations;
    const double t = std::log(static_cast<double>(iterations) / sliderMinIterations) / std::log(ratio);
    return static_cast<int>(std::lround(t * sliderSteps));
}

}  // namespace

QtMainWindow::QtMainWindow(QWidget* parent) : QMainWindow(parent), m_centralWidget(new QMandelbrotWidget(this))
{
    ui.setupUi(this);
    setWindowTitle("Mandelbrot (Qt6)");
    CreateIterationsPanel();
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
    precisionGroup->addAction(ui.actionPerturbation);
    precisionGroup->setExclusive(true);
    connect(ui.actionPrecisionAuto, &QAction::triggered, this, &QtMainWindow::onActionPrecision);
    connect(ui.actionFixedPoint128, &QAction::triggered, this, &QtMainWindow::onActionPrecision);
    connect(ui.actionDouble, &QAction::triggered, this, &QtMainWindow::onActionPrecision);
    connect(ui.actionPerturbation, &QAction::triggered, this, &QtMainWindow::onActionPrecision);

    // Palette actions
    QActionGroup* paletteGroup = new QActionGroup(this);
    paletteGroup->addAction(ui.actionPaletteGrey);
    paletteGroup->addAction(ui.actionPaletteGradient);
    paletteGroup->addAction(ui.actionPaletteHistorgram);
    paletteGroup->addAction(ui.actionPaletteVivid);
    paletteGroup->setExclusive(true);
    connect(ui.actionPaletteGrey, &QAction::triggered, [this]() { m_centralWidget->setPaletteType(QMandelbrotWidget::palGrey); });
    connect(ui.actionPaletteGradient, &QAction::triggered, [this]() { m_centralWidget->setPaletteType(QMandelbrotWidget::palGradient); });
    connect(ui.actionPaletteHistorgram, &QAction::triggered, [this]() { m_centralWidget->setPaletteType(QMandelbrotWidget::palHistogram); });
    connect(ui.actionPaletteVivid, &QAction::triggered, [this]() { m_centralWidget->setPaletteType(QMandelbrotWidget::palVivid); });
    
    // smooth transitions toggle
    connect(ui.actionSmoothTranisitions, &QAction::toggled, m_centralWidget, &QMandelbrotWidget::setSmoothTransitions);

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
    
    // default palette
    ui.actionPaletteGradient->setChecked(true);

    // default smooth coloring
    ui.actionSmoothTranisitions->setChecked(true);
}

/**
 * @brief Build the iterations slider panel and lay it out left of the fractal widget.
 *
 * The slider spans [0, sliderSteps] and maps to [sliderMinIterations, sliderMaxIterations]
 * on a log scale. Tracking is off, so a drag renders once on release while the label follows
 * the handle; wheel steps and groove clicks still apply immediately.
 */
void QtMainWindow::CreateIterationsPanel()
{
    _iterSlider = new QSlider(Qt::Vertical);
    _iterSlider->setRange(0, sliderSteps);
    _iterSlider->setSingleStep(sliderSingleStep);
    _iterSlider->setPageStep(sliderPageStep);
    _iterSlider->setTracking(false);
    // the main window handles the arrow and +/- keys for navigation, so the slider must not take focus
    _iterSlider->setFocusPolicy(Qt::NoFocus);
    _iterSlider->setToolTip("Iteration limit (logarithmic scale)");
    _iterSlider->setValue(IterationsToSliderPos(QMandelbrotWidget::min_iterations));

    _iterLabel = new QLabel;
    _iterLabel->setAlignment(Qt::AlignCenter);
    _iterLabel->setToolTip("Approximate iteration limit; the status bar shows the exact value");
    // size for the widest value so the panel, and with it the fractal, never resizes as the text changes
    _iterLabel->setFixedWidth(_iterLabel->fontMetrics().horizontalAdvance(QString::number(sliderMaxIterations)) + 8);
    UpdateIterationsLabel(_iterSlider->value());

    QVBoxLayout* panelLayout = new QVBoxLayout;
    panelLayout->setContentsMargins(4, 4, 4, 4);
    panelLayout->addWidget(_iterLabel);
    panelLayout->addWidget(_iterSlider, 1, Qt::AlignHCenter);

    QWidget* container = new QWidget(this);
    QHBoxLayout* containerLayout = new QHBoxLayout(container);
    containerLayout->setContentsMargins(0, 0, 0, 0);
    containerLayout->setSpacing(0);
    containerLayout->addLayout(panelLayout);
    containerLayout->addWidget(m_centralWidget, 1);
    setCentralWidget(container);

    connect(_iterSlider, &QSlider::sliderMoved, this, &QtMainWindow::UpdateIterationsLabel);
    connect(_iterSlider, &QSlider::valueChanged, this, &QtMainWindow::OnIterationsSliderChanged);
}

/**
 * @brief Apply an iteration limit picked with the slider.
 *
 * Only user input reaches here; programmatic moves go through SyncIterationsSlider(),
 * which blocks signals.
 *
 * @param position Slider position in [0, sliderSteps].
 */
void QtMainWindow::OnIterationsSliderChanged(int position)
{
    UpdateIterationsLabel(position);

    // the slider is in control now, so no menu item (Auto included) stays checked
    if (QAction* act = iterGroup->checkedAction()) {
        act->setChecked(false);
    }

    m_centralWidget->setMaximumIterations(SliderPosToIterations(position));
}

/**
 * @brief Move the slider to the position nearest an iteration limit without applying it.
 * @param iterations Iteration limit to show.
 */
void QtMainWindow::SyncIterationsSlider(int64_t iterations)
{
    // don't pull the handle away from the user mid-drag
    if (_iterSlider->isSliderDown()) {
        return;
    }

    const QSignalBlocker blocker(_iterSlider);
    _iterSlider->setValue(IterationsToSliderPos(iterations));
    UpdateIterationsLabel(_iterSlider->value());
}

/**
 * @brief Show the iteration limit at a slider position in the slider label.
 * @param position Slider position in [0, sliderSteps].
 */
void QtMainWindow::UpdateIterationsLabel(int position)
{
    _iterLabel->setText(QString::number(SliderPosToIterations(position)));
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
    else if (act == ui.actionPerturbation)
        m_centralWidget->setPrecision(QMandelbrotWidget::Precision::Perturbation);
}

/**
 * @brief Parse the selected iteration action and forward to the widget.
 *
 * If the action text parses as an integer, that value is used directly and the
 * slider moves to match. Otherwise, automatic iteration scaling is enabled
 * (maxIter = 0) and the slider follows the computed limit from onRenderDone().
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
        SyncIterationsSlider(iter);
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
 * (e.g. "2^42"). Otherwise, a decimal value is shown. While Auto iterations is
 * active, the slider also moves to the limit used for this frame.
 *
 * @param stats Frame statistics containing render time, zoom, size, and iterations.
 */
void QtMainWindow::onRenderDone(FrameStats stats)
{
    // in Auto mode the slider follows the limit derived from the zoom level
    if (ui.actionIterAuto->isChecked()) {
        SyncIterationsSlider(stats.max_iterations);
    }

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
