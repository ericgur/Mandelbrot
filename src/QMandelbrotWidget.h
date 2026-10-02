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

#include <QWidget>
#include <QChronoTimer>

// include fixed point implementation from project root
#include "fp128/fixed_point128.h"

/**
 * @brief Integer bits of fp128_t, not counting the sign.
 *
 * Every bit taken from the integer part goes to the fraction and buys one more zoom level, so
 * this is as small as the values the renderer holds allow: 3 bits hold [-8, 8). The iterate
 * after a step stays within 4 + 2.83 of the origin and pixel coordinates within 7, both of which
 * fit. The squared modulus of the iterate that escapes, up to 46.6, does not; the escape test in
 * QMandelbrotWidget.cpp is written to stay exact when it wraps. Any value in [3, 63] renders
 * correctly, each one level shallower than the one below it.
 */
inline constexpr int32_t fp128IntBits = 3;

/// @brief 128-bit fixed-point type with fp128IntBits integer bits, a sign bit and 127 - fp128IntBits fraction bits.
typedef fp128::fixed_point128<fp128IntBits> fp128_t;

/**
 * @struct Complex128
 * @brief A complex number with fp128 parts, such as the Julia constant.
 *
 * std::complex is only specified for the built-in floating point types, so the parts are
 * two plain fp128_t members.
 */
struct Complex128 {
    fp128_t real {};  ///< Real part.
    fp128_t imag {};  ///< Imaginary part.

    /**
     * @brief Compare both parts exactly.
     * @param rhs Value to compare with.
     * @return True if both parts are equal.
     */
    bool operator==(const Complex128& rhs) const = default;
};

/**
 * @struct FrameStats
 * @brief Statistics collected after each frame render.
 */
struct FrameStats {
    uint32_t render_time_ms {};  ///< Wall-clock render time in milliseconds.
    int32_t log2Zoom {};         ///< Zoom level as a power of 2.
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
 * for deeper zooms up to 2^logMaxZoom), and provides four color palette modes
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
    /// Log2 of maximum zoom: the deepest level at which neighboring pixels of a 3840 pixel wide image
    /// are still at least one fp128_t LSB apart (114 with 3 integer bits).
    static inline constexpr int32_t logMaxZoom = fp128_t::F - 10;
    static inline constexpr int32_t logMinZoom = 0;         ///< Log2 of minimum zoom (x1).
    /// Largest magnitude of either part of the view center. The Mandelbrot set lies within 2 of
    /// the origin, and so does every Julia set that is more than dust, so nothing past it is lost.
    static inline constexpr int32_t maxCenterMagnitude = 2;
    /// setMaximumIterations() value that turns on Auto, which scales the iteration limit with the zoom.
    static inline constexpr int64_t auto_iterations = 0;
    /// Lowest fixed iteration limit the UI offers; max_iterations is the highest.
    static inline constexpr int64_t min_fixed_iterations = 64;

    /** @brief Rendering precision modes. */
    enum class Precision {
        Auto,           ///< Double up to zoom 2^44, then FixedPoint128.
        Double,         ///< Force IEEE 754 double precision.
        FixedPoint128,  ///< Force 128-bit fixed-point precision.
        Perturbation    ///< Perturbation theory with fp128 reference + double deltas (Mandelbrot only).
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
     * @brief Render one frame into the internal image cache without painting to the screen.
     *
     * Runs the exact pipeline paintEvent() drives - buffer (re)allocation, the escape-time
     * calculation at the active precision, palette/histogram rebuild and the iteration-to-RGB
     * conversion - and stops short of blitting the result to the widget surface. The widget
     * does not need to be visible, which is what lets the benchmark mode drive it headless.
     *
     * The frame is rendered at the widget's current size(), so callers should resize() first.
     *
     * @return Wall-clock render time in milliseconds, with sub-millisecond resolution.
     */
    [[nodiscard]] double renderOffscreen();

    /**
     * @brief Point the view at a specific location in the complex plane at a given zoom.
     *
     * The default view spans x in [-2.5, 2.5], and each zoom step halves that span, so the
     * horizontal half-width becomes 2.5 / 2^log2Zoom. The vertical extent follows from the
     * widget aspect ratio about @p centerY.
     *
     * @param centerX Real part of the view center, clamped to [-maxCenterMagnitude, maxCenterMagnitude].
     * @param centerY Imaginary part of the view center, clamped the same way.
     * @param log2Zoom Log2 of the zoom level, clamped to [logMinZoom, logMaxZoom].
     */
    void setView(const fp128_t& centerX, const fp128_t& centerY, int32_t log2Zoom);

    /**
     * @brief Get the current Julia set constant.
     * @return The complex constant C used for Julia set rendering.
     */
    [[nodiscard]] Complex128 juliaConstant() const { return _juliaConstant; }

    /**
     * @brief Get the Julia constant used until another one is set.
     * @return 0.285 + 0.01i, exact to fp128 precision.
     */
    [[nodiscard]] static Complex128 defaultJuliaConstant();

    /**
     * @brief Get the real part of the view center.
     * @return The view center's real part, at full fp128 precision.
     */
    [[nodiscard]] fp128_t viewCenterX() const { return _centerX; }

    /**
     * @brief Get the imaginary part of the view center.
     * @return The view center's imaginary part, at full fp128 precision.
     */
    [[nodiscard]] fp128_t viewCenterY() const { return _centerY; }

    /**
     * @brief Get the zoom level as a power of 2, in the form setView() takes it.
     * @return Log2 of the zoom level, in [logMinZoom, logMaxZoom].
     */
    [[nodiscard]] int32_t log2Zoom() const { return _logZoomLevel; }

    /**
     * @brief Get the fractal set being rendered.
     * @return The active set type.
     */
    [[nodiscard]] set_type_t setType() const { return _setType; }

    /**
     * @brief Get the iteration limit setting, in the form setMaximumIterations() takes it.
     * @return The fixed iteration limit, or auto_iterations while Auto is on.
     */
    [[nodiscard]] int64_t maximumIterations() const { return _autoIterations ? auto_iterations : _maxIter; }

    /**
     * @brief Get the iteration limit Auto picks at a zoom level.
     *
     * Scales linearly with log2 of the zoom, from min_iterations at 1x to max_iterations at 2^logMaxZoom.
     *
     * @param log2Zoom Log2 of the zoom level; values below 0 count as 0.
     * @return The iteration limit.
     */
    [[nodiscard]] static int64_t autoIterationLimit(int32_t log2Zoom);

    /**
     * @brief Switch between Mandelbrot and Julia set rendering.
     * @param type The fractal set type to render.
     */
    void setSetType(set_type_t type);

    /**
     * @brief Query whether OpenMP parallelization is enabled.
     * @return True if OpenMP is active.
     */
    bool openMp() const { return _useOpenMP; }

    /**
     * @brief Mark the fractal for recomputation on the next frame, and the color table too if asked.
     *
     * Passing false leaves the color table as it is; it never marks a stale table valid. A table
     * left stale by an iteration limit or palette change must still be rebuilt, or the frame
     * reads it past its end with the new limit.
     *
     * @param invalidateColorTable True to rebuild the color table as well.
     */
    virtual void invalidate(bool invalidateColorTable = true)
    {
        setFractalDataValid(false);
        if (invalidateColorTable) {
            setColorTableValid(false);
        }
        update();
    }

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
    void setJuliaConstant(const Complex128& c);

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
    void setOpenMp(bool enable) { _useOpenMP = enable; }

    /**
     * @brief Modify the palette used to draw the fractal.
     * @param palette The palette type to use.
     */
    void setPaletteType(palette_t palette);

    /**
     * @brief Enable or disable smooth transitions between palette colors.
     * @param enable True to enable smooth transitions of palette entries, false to disable.
     */
    void setSmoothTransitions(bool enable);

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
     * max_iterations (at zoom 2^logMaxZoom).
     *
     * @return The computed iteration limit.
     */
    int64_t calcAutoIterationLimits();

    /**
     * @brief Render one frame into _imageCache.
     *
     * The shared body of paintEvent() and renderOffscreen(): (re)allocates the image and
     * iteration buffers when the widget size changed, refreshes the auto-iteration limit,
     * runs the escape-time calculation at the active precision, rebuilds the palette when
     * it went stale and converts the iteration buffer to RGB. Nothing here touches the
     * widget surface, so it is safe to call on a widget that was never shown.
     *
     * @return Wall-clock render time in nanoseconds.
     */
    [[nodiscard]] int64_t RenderFrame();

    inline bool fractalDataValid() const { return _fractalDataValid; }
    inline void setFractalDataValid(bool valid = true) { _fractalDataValid = valid; }
    inline bool colorTableValid() const { return _colorTableValid; }
    inline void setColorTableValid(bool valid = true) { _colorTableValid = valid; }

    // View state
    fp128_t _centerX, _centerY;             ///< View center in the complex plane, each part in [-maxCenterMagnitude, maxCenterMagnitude].
    int32_t _logZoomLevel = 0;              ///< Zoom level as a power of 2, in [logMinZoom, logMaxZoom]; 0 = default.
    int64_t _maxIter = 128;                 ///< Current maximum iteration count.
    bool _autoIterations = false;           ///< True if iterations scale with zoom.

    // Image and buffers
    QImage _imageCache;                    ///< Cached rendered QImage.
    float* _iterations = nullptr;          ///< Per-pixel iteration count buffer.
    bool _fractalDataValid = false;        ///< False when the fractal needs recomputation.
    bool _colorTableValid = false;         ///< False when color tables need to rebuild, e.g. when changing from palette to another.

    // Color / palette data
    QVector<QRgb> _colorTable;             ///< Color lookup table indexed by iteration count.
    int* _histogram = nullptr;             ///< Iteration frequency histogram for palette equalization.
    float _hsvOffset = 0;                  ///< HSV hue rotation offset for animation.
    palette_t _paletteType = palGradient;  ///< Active color palette type.
    bool _smoothLevel = true;              ///< Enable smooth (fractional) iteration coloring.

    // UI flags
    Precision _precision = Precision::Auto;  ///< Active rendering precision mode.
    bool _animate = false;                   ///< True if palette animation is running.

    // Set type and Julia constants
    set_type_t _setType = stMandelbrot;                  ///< Active fractal set type.
    Complex128 _juliaConstant = defaultJuliaConstant();  ///< Julia set complex constant.

    // Timer for animation
    QChronoTimer _timer;  ///< Timer driving palette animation ticks.

    // OpenMP support
    bool _useOpenMP = true;  ///< OpenMP parallelization toggle.

    // Helpers

    /** @brief Initialize view bounds to the default complex plane region. */
    void SetDefaultValues();

    /**
     * @brief Move the view center, keeping each part within maxCenterMagnitude.
     * @param centerX Real part of the view center.
     * @param centerY Imaginary part of the view center.
     */
    void SetViewCenter(const fp128_t& centerX, const fp128_t& centerY);

    /**
     * @brief Get half the width of the view in the complex plane.
     * @return 2.5 / 2^_logZoomLevel, exact at every zoom level.
     */
    [[nodiscard]] fp128_t ViewHalfWidth() const;

    /**
     * @brief Run the escape-time calculation for the current view at the active precision.
     *
     * The view keeps its center and horizontal span at any image size; the vertical span
     * follows from the image's aspect ratio.
     *
     * @param pIterations Output buffer for per-pixel iteration counts.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     */
    void CalcIterations(float* pIterations, int64_t width, int64_t height);

    /**
     * @brief Zoom by a power of 2 and center the view on a screen point.
     *
     * The new zoom level is clamped to [logMinZoom, logMaxZoom]; when the clamp leaves it
     * unchanged, the view is left as it is.
     *
     * @param point Screen coordinates of the new view center.
     * @param logZoomDelta Zoom steps of 2x each (positive to zoom in, negative to zoom out).
     */
    void OnZoomChange(const QPoint& point, int32_t logZoomDelta);

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
     *
     * Dispatches to the templated implementation based on the active set type
     * so the inner iteration loop is fully specialized (no per-pixel branching
     * between Mandelbrot and Julia code paths).
     *
     * @param pIterations Output buffer for per-pixel iteration counts.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     * @param x0 Left edge of the view in the complex plane.
     * @param dx Horizontal step per pixel.
     * @param y0 Top edge of the view in the complex plane.
     * @param dy Vertical step per pixel.
     */
    void CalcIterationsDouble(float* pIterations, int64_t width, int64_t height, double x0, double dx, double y0, double dy);

    /**
     * @brief Render the fractal using 128-bit fixed-point precision.
     *
     * Dispatches to the templated implementation based on the active set type.
     *
     * @param pIterations Output buffer for per-pixel iteration counts.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     * @param centerX Real part of the view center.
     * @param centerY Imaginary part of the view center.
     * @param halfWidth Half the view width in the complex plane.
     */
    void CalcIterationsFP128(float* pIterations, int64_t width, int64_t height, const fp128_t& centerX, const fp128_t& centerY, const fp128_t& halfWidth);

    /**
     * @brief Templated double-precision render specialized on set type.
     * @tparam IsJulia True for Julia, false for Mandelbrot.
     */
    template<bool IsJulia>
    void CalcIterationsDoubleImpl(float* pIterations, int64_t width, int64_t height, double x0, double dx, double y0, double dy);

    /**
     * @brief Templated 128-bit fixed-point render specialized on set type.
     * @tparam IsJulia True for Julia, false for Mandelbrot.
     */
    template<bool IsJulia>
    void CalcIterationsFP128Impl(float* pIterations, int64_t width, int64_t height, const fp128_t& centerX, const fp128_t& centerY, const fp128_t& halfWidth);

    /**
     * @brief Render the Mandelbrot set using perturbation theory.
     *
     * Computes one high-precision reference orbit (fp128) at the view center,
     * then iterates each pixel as a small @c double delta around that orbit.
     * Pixels that glitch (Pauldelbrot's criterion: |Z+Δ|² ≪ |Z|²) or that
     * iterate past the reference's escape point fall back to per-pixel fp128.
     *
     * Falls back to CalcIterationsFP128Impl for Julia (perturbation needs a
     * different reformulation that isn't implemented here).
     *
     * @param pIterations Output buffer for per-pixel iteration counts.
     * @param width Image width in pixels.
     * @param height Image height in pixels.
     * @param centerX Real part of the view center.
     * @param centerY Imaginary part of the view center.
     * @param halfWidth Half the view width in the complex plane.
     */
    void CalcIterationsPerturbation(float* pIterations, int64_t width, int64_t height, const fp128_t& centerX, const fp128_t& centerY,
                                    const fp128_t& halfWidth);

    /**
     * @brief Compute a single Mandelbrot pixel at full fp128 precision.
     *
     * Used by the perturbation path as a fallback for glitched pixels and
     * pixels whose orbit outlives the reference. Inlined helper rather than
     * a separate template to avoid touching the main fp128 dispatcher.
     *
     * @param cx Pixel real coordinate in the complex plane.
     * @param cy Pixel imaginary coordinate in the complex plane.
     * @return Smooth or integer iteration count, depending on _smoothLevel.
     */
    float CalcSinglePixelFP128(fp128_t cx, fp128_t cy);
};
