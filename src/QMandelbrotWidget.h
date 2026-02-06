/**
 * @file QMandelbrotWidget.h
 * @brief Core fractal rendering engine widget for the qMandelbrot application.
 *
 * Implements the escape-time algorithm for both Mandelbrot and Julia sets
 * with dual-precision rendering (IEEE 754 doubles and custom 128-bit
 * fixed-point arithmetic). Supports multiple color palettes with animation,
 * smooth iteration interpolation, and OpenMP-parallelized scanline rendering.
 */

#pragma once

#include <complex>
#include <QWidget>
#include <QChronoTimer>

// include fixed point implementation from project root
#include "fixed_point128.h"
using namespace fp128;

/// @brief 128-bit fixed-point type with 8 integer bits and 120 fractional bits.
typedef fixed_point128<8> fp128_t;

/**
 * @struct FrameStats
 * @brief Statistics collected after each frame render.
 */
struct FrameStats {
    uint32_t render_time_ms {};  ///< Wall-clock render time in milliseconds.
    float zoom {};               ///< Current zoom level multiplier.
    QSize size {};               ///< Rendered image dimensions in pixels.
    int32_t max_iterations {};   ///< Maximum iteration count used for this frame.
};

/**
 * @class QMandelbrotWidget
 * @brief Widget that renders Mandelbrot and Julia set fractals.
 *
 * This is the core rendering engine of the application. It computes
 * the escape-time algorithm per-pixel, supports dual-precision rendering
 * (double-precision IEEE 754 up to zoom 2^44, then 128-bit fixed-point
 * for deeper zooms up to 2^113), and provides four color palette modes
 * (grey, gradient, vivid, histogram-equalized) with optional animation.
 *
 * Rendering is parallelized across scanlines using OpenMP with dynamic
 * scheduling. Smooth coloring is achieved via logarithmic interpolation
 * of fractional iteration counts.
 *
 * @par Mouse Controls
 * - Left click: Zoom in 2x (4x with Ctrl, 8x with Ctrl+Shift)
 * - Right click: Zoom out 2x (4x with Ctrl, 8x with Ctrl+Shift)
 * - Middle click: Reset to default view
 */
class QMandelbrotWidget : public QWidget
{
    Q_OBJECT
    Q_DISABLE_COPY_MOVE(QMandelbrotWidget)
public:
    static inline constexpr int64_t max_iterations = 2500;  ///< Upper bound for iteration count.
    static inline constexpr int64_t min_iterations = 128;   ///< Lower bound for iteration count.
    static inline constexpr double logMaxZoom = 113.0;      ///< Log2 of maximum zoom (128-bit fixed-point limit).

    /** @brief Rendering precision modes. */
    enum class Precision {
        Auto,          ///< Double up to zoom 2^44, then FixedPoint128.
        Double,        ///< Force IEEE 754 double precision.
        FixedPoint128  ///< Force 128-bit fixed-point precision.
    };

    /** @brief Color palette types. */
    enum palette_t {
        palGrey,      ///< Smooth greyscale gradient.
        palGradient,  ///< Progressive RGB channel shifts.
        palVivid,     ///< 6-segment HSV rainbow cycle.
        palHistogram  ///< Histogram-equalized HSV distribution.
    };

    /** @brief Fractal set types. */
    enum set_type_t {
        stMandelbrot,  ///< Standard Mandelbrot set.
        stJulia,       ///< Julia set with configurable constant.
        stCount        ///< Sentinel value for set type count.
    };

    /**
     * @brief Construct the rendering widget.
     * @param parent Optional parent widget.
     */
    explicit QMandelbrotWidget(QWidget* parent = nullptr);

    /** @brief Destructor. Frees iteration and histogram buffers. */
    virtual ~QMandelbrotWidget();

    /**
     * @brief Export the current view as a PNG image at the specified resolution.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     */
    void saveImage(int width, int height);

    /**
     * @brief Get the current Julia set constant.
     * @return The complex constant C used for Julia set rendering.
     */
    std::complex<double> juliaConstant() const { return m_JuliaConstant; }

    /**
     * @brief Switch between Mandelbrot and Julia set rendering.
     * @param type The fractal set type to render.
     */
    void setSetType(set_type_t type);

    /**
     * @brief Query whether OpenMP parallelization is enabled.
     * @return True if OpenMP is active.
     */
    bool openMp() const { return m_UseOpenMP; }

signals:
    /**
     * @brief Emitted after each frame render with timing and view statistics.
     * @param stats The frame statistics.
     */
    void renderDone(FrameStats stats);

public slots:
    /**
     * @brief Set the Julia set complex constant and trigger a re-render.
     * @param c The new complex constant value.
     */
    void setJuliaConstant(const std::complex<double>& c);

    /** @brief Reset the view to the default bounds and zoom level. */
    void resetView();

    /** @brief Zoom in 2x from the center of the viewport. */
    void zoomIn();

    /** @brief Zoom out 2x from the center of the viewport. */
    void zoomOut();

    /** @brief Advance palette animation by one tick. */
    void animationTick();

    /**
     * @brief Pan the view horizontally.
     * @param amount Fraction of the viewport width to pan (positive = right).
     */
    void panHorizontal(double amount);

    /**
     * @brief Pan the view vertically.
     * @param amount Fraction of the viewport height to pan (positive = down).
     */
    void panVertical(double amount);

    /**
     * @brief Enable or disable palette color cycling animation.
     * @param animate True to start animation, false to stop and reset.
     */
    void setAnimatePalette(bool animate);

    /**
     * @brief Set the rendering precision mode.
     * @param p The precision mode to use.
     */
    void setPrecision(Precision p);

    /**
     * @brief Set the maximum iteration count.
     * @param maxIter Iteration limit, or 0 to enable automatic scaling.
     */
    void setMaximumIterations(int64_t maxIter);

    /**
     * @brief Enable or disable OpenMP parallel rendering.
     * @param enable True to enable multi-threaded rendering.
     */
    void setOpenMp(bool enable) { m_UseOpenMP = enable; }

protected:
    /**
     * @brief Render the fractal image and paint it to the widget.
     * @param event Paint event (unused).
     */
    void paintEvent(QPaintEvent* event) override;

    /**
     * @brief Handle window resize by marking the view for recomputation.
     * @param event Resize event.
     */
    void resizeEvent(QResizeEvent* event) override;

    /**
     * @brief Handle mouse clicks for zooming and view reset.
     * @param event Mouse event with button and modifier information.
     */
    void mousePressEvent(QMouseEvent* event) override;

private:
    /**
     * @brief Compute iteration limits based on the current zoom level.
     *
     * Linearly interpolates between min_iterations (at zoom 1) and
     * max_iterations (at zoom 2^113).
     *
     * @return The computed iteration limit.
     */
    int64_t calcAutoIterationLimits();

    // View state
    fp128_t m_Xmin, m_Xmax, m_Ymin, m_Ymax;  ///< View bounds in the complex plane.
    double m_ZoomLevel;                      ///< Current zoom multiplier (1.0 = default).
    double m_ZoomIncrement = 2.0;            ///< Zoom step factor per click.
    int64_t m_MaxIter = 128;                 ///< Current maximum iteration count.
    bool m_AutoIterations = false;           ///< True if iterations scale with zoom.

    // Image and buffers
    QImage m_ImageCache;            ///< Cached rendered QImage.
    float* m_Iterations = nullptr;  ///< Per-pixel iteration count buffer.
    bool m_NeedToRecompute = true;  ///< True when the fractal needs recomputation.

    // Color / palette data
    QVector<QRgb> m_ColorTable;             ///< Color lookup table indexed by iteration count.
    int* m_Histogram = nullptr;             ///< Iteration frequency histogram for palette equalization.
    float m_HsvOffset = 0;                  ///< HSV hue rotation offset for animation.
    palette_t m_PaletteType = palGradient;  ///< Active color palette type.
    bool m_SmoothLevel = true;              ///< Enable smooth (fractional) iteration coloring.

    // UI flags
    Precision m_Precision = Precision::Auto;  ///< Active rendering precision mode.
    bool m_Animate = false;                   ///< True if palette animation is running.

    // Set type and Julia constants
    set_type_t m_SetType = stMandelbrot;                 ///< Active fractal set type.
    std::complex<double> m_JuliaConstant {0.285, 0.01};  ///< Julia set complex constant.

    // Timer for animation
    QChronoTimer m_Timer;  ///< Timer driving palette animation ticks.

    // OpenMP support
    bool m_UseOpenMP = true;  ///< OpenMP parallelization toggle.

    // Helpers

    /** @brief Initialize view bounds to the default complex plane region. */
    void SetDefaultValues();

    /** @brief Adjust Y bounds to maintain the correct aspect ratio. */
    void SetAspectRatio();

    /**
     * @brief Zoom the view centered on a screen point.
     * @param point Screen coordinates of the zoom center.
     * @param zoomMultiplier Zoom factor (>1 to zoom in, <1 to zoom out).
     */
    void OnZoomChange(const QPoint& point, double zoomMultiplier);

    // Rendering helpers

    /** @brief Build the color lookup table for the current palette type. */
    void CreateColorTables();

    /**
     * @brief Generate a histogram-equalized HSV color table.
     * @param offset HSV hue rotation offset for animation.
     */
    void CreateColorTableFromHistogram(float offset);

    /**
     * @brief Build an iteration frequency histogram from the iteration buffer.
     * @param pIterations Per-pixel iteration count buffer.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     */
    void CreateHistogram(const float* pIterations, int64_t width, int64_t height);

    /**
     * @brief Convert the iteration buffer to an RGB QImage using the color table.
     * @param img Output QImage (must be Format_RGB32).
     * @param pIterations Per-pixel iteration count buffer.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     */
    void CreateDibFromIterations(QImage& img, const float* pIterations, int64_t width, int64_t height);

    /**
     * @brief Render the fractal using IEEE 754 double precision.
     * @param pIterations Output buffer for per-pixel iteration counts.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     * @param x0 Left edge of the view in the complex plane.
     * @param dx Horizontal step per pixel.
     * @param y0 Top edge of the view in the complex plane.
     * @param dy Vertical step per pixel.
     */
    void DrawImageDouble(float* pIterations, int64_t width, int64_t height, double x0, double dx, double y0, double dy);

    /**
     * @brief Render the fractal using 128-bit fixed-point precision.
     * @param pIterations Output buffer for per-pixel iteration counts.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     * @param x0 Left edge of the view in the complex plane.
     * @param dx Horizontal step per pixel.
     * @param y0 Top edge of the view in the complex plane.
     * @param dy Vertical step per pixel.
     */
    void DrawImageFixedPoint128(float* pIterations, int64_t width, int64_t height, fp128_t x0, fp128_t dx, fp128_t y0, fp128_t dy);
};
