/**
 * @file QMandelbrotWidget.cpp
 * @brief Implementation of the QMandelbrotWidget fractal rendering engine.
 *
 * Contains the escape-time rendering loops (double and 128-bit fixed-point),
 * color palette generation (grey, gradient, vivid, histogram-equalized),
 * smooth iteration interpolation, view management, and mouse/keyboard input handling.
 */

#include "pch.h"
#include "QMandelbrotWidget.h"
#include <QPainter>
#include <QFileDialog>
#include <QMouseEvent>
#include <QElapsedTimer>
#include <cmath>

using namespace std::chrono_literals;

constexpr uint64_t MAX_DOUBLE_ZOOM_LEVEL = 1ull << 44;

/**
 * @brief Blend two QRgb colors using alpha interpolation.
 *
 * Performs fast integer-based alpha blending without floating-point arithmetic.
 * Alpha of 0 returns c1, alpha of 256 returns c2.
 *
 * @param c1 First color (alpha = 0).
 * @param c2 Second color (alpha = 256).
 * @param alpha Blending factor in [0, 256].
 * @return The blended QRgb color.
 */
static inline QRgb blendAlphaQRgb(QRgb c1, QRgb c2, uint32_t alpha)
{
    uint32_t rb1 = (0x100 - alpha) * ((c1 & 0xFF00FF));
    uint32_t rb2 = alpha * ((c2 & 0xFF00FF));
    uint32_t g1 = (0x100 - alpha) * ((c1 & 0x00FF00));
    uint32_t g2 = alpha * ((c2 & 0x00FF00));
    uint32_t rb = (((rb1 + rb2) >> 8) & 0xFF00FF);
    uint32_t g = (((g1 + g2) >> 8) & 0x00FF00);
    return (QRgb)(rb | g);
}

QMandelbrotWidget::QMandelbrotWidget(QWidget* parent) : QWidget(parent), _timer(this)
{
    SetDefaultValues();
    CreateColorTables();
    connect(&_timer, &QChronoTimer::timeout, this, &QMandelbrotWidget::animationTick);
}

QMandelbrotWidget::~QMandelbrotWidget()
{
    delete[] _iterations;
    delete[] _histogram;
}

/**
 * @brief Initialize view bounds to the default complex plane region.
 *
 * Sets the X range to [-2.5, 2.5] and zoom to 1x. If the INITIAL_POINT
 * macro is defined, uses preset coordinates for a specific deep-zoom location.
 */
void QMandelbrotWidget::SetDefaultValues()
{
#ifdef INITIAL_POINT
    _xmin = x_init[0];
    _xmax = x_init[1];
    _ymin = y_init[0];
    _ymax = y_init[1];
    _zoomLevel = zoom_init;
#else
    _xmax = 2.5;
    _xmin = -_xmax;
    _ymax = _ymin = 0;
    _zoomLevel = 1;
#endif
    SetAspectRatio();
}

/**
 * @brief Adjust Y bounds to preserve the correct aspect ratio.
 *
 * Recomputes _ymin and _ymax based on the current X bounds and
 * the widget's width/height ratio, keeping the Y center unchanged.
 */
void QMandelbrotWidget::SetAspectRatio()
{
    QSize s = size();

    // use _xmin, _xmax and m_rect to determine _ymin and _ymax
    // check if window created
    if (0 == s.height() || 0 == s.width())
        return;

    fp128_t ratio = (double)s.height() / s.width();
    fp128_t ysize = (_xmax - _xmin) * (ratio >> 1);
    _ymin = ((_ymax + _ymin) >> 1) - ysize;
    _ymax = _ymin + (ysize << 1);
}

/**
 * @brief Zoom the view centered on a specific screen coordinate.
 *
 * Preserves the complex-plane coordinate under the given pixel while
 * scaling the view bounds by the specified multiplier.
 *
 * @param point Screen coordinates of the zoom center.
 * @param zoomMultiplier Zoom factor (>1 to zoom in, <1 to zoom out).
 */
void QMandelbrotWidget::OnZoomChange(const QPoint& point, double zoomMultiplier)
{
    QSize s = size();
    static fp128_t one = 1;

    // fix y coords
    fp128_t alpha = (double)(point.y()) / (s.height() - 1);
    fp128_t quarter = (_ymax - _ymin) * (1.0 / (zoomMultiplier * 2.0));
    fp128_t center = alpha * _ymax + (one - alpha) * _ymin;

    _ymin = center - quarter;
    _ymax = center + quarter;

    // fix x coords
    alpha = (double)(point.x()) / (s.width() - 1);
    quarter = (_xmax - _xmin) * (1.0 / (zoomMultiplier * 2.0));
    center = alpha * _xmax + (one - alpha) * _xmin;
    _xmin = center - quarter;
    _xmax = center + quarter;

    invalidate(false);
}

/**
 * @brief Build the color lookup table for the current palette type.
 *
 * Generates _maxIter + 1 color entries based on the active palette:
 * - @b Grey: Smooth greyscale from white to dark grey.
 * - @b Vivid: 6-segment HSV rainbow cycle.
 * - @b Gradient: Progressive RGB channel shifts (+3, +5, -3).
 *
 * Entry 0 is always white (escaped at iteration 0) and entry _maxIter
 * is always black (inside the set).
 */
void QMandelbrotWidget::CreateColorTables()
{
    _colorTable.clear();
    _colorTable.resize(_maxIter + 1);

    if (_paletteType == palGrey) {
        for (int64_t i = 1; i <= _maxIter; ++i) {
            int c = 255 - (int)(215.0f * (float)i / (float)_maxIter);
            _colorTable[i] = qRgb(c, c, c);
        }
    } else if (_paletteType == palVivid) {
        float step = 13.0f / 256.f;
        for (int i = 1; i < _maxIter; ++i) {
            float h = step * i;
            float x = (1.0f - fabs(fmodf(h, 2) - 1.0f));
            switch ((int)floorf(h) % 6) {
            case 0:  // 0-60
                _colorTable[i] = qRgb(255, int(x * 255), 0);
                break;
            case 1:  // 60-120
                _colorTable[i] = qRgb(int(x * 255), 255, 0);
                break;
            case 2:  // 120-180
                _colorTable[i] = qRgb(0, 255, int(x * 255));
                break;
            case 3:  // 180-240
                _colorTable[i] = qRgb(0, int(x * 255), 255);
                break;
            case 4:  // 240-300
                _colorTable[i] = qRgb(int(x * 255), 0, 255);
                break;
            case 5:  // 300-360
                _colorTable[i] = qRgb(255, 0, int(x * 255));
                break;
            }
        }
    } else if (_paletteType == palGradient) {
        uint32_t r = 0;
        uint32_t g = 20;
        uint32_t b = 255;

        for (int64_t i = 1; i <= _maxIter; ++i) {
            r = (r + 3) & 0xFF;
            g = (g + 5) & 0xFF;
            b = (b - 3) & 0xFF;
            _colorTable[i] = qRgb(r, g, b);
        }
    }

    _colorTable[0] = qRgb(255, 255, 255);
    _colorTable[_maxIter] = qRgb(0, 0, 0);
    setColorTableValid();
}

/**
 * @brief Generate a histogram-equalized HSV color table.
 *
 * Distributes hues across iteration counts proportionally to their frequency
 * in the histogram, clamped by a threshold to prevent a single iteration
 * from dominating the palette. The offset parameter rotates the hue wheel
 * for animation.
 *
 * @param offset HSV hue rotation offset in [0, 1).
 */
void QMandelbrotWidget::CreateColorTableFromHistogram(float offset)
{
    // Initialize to the -1.0 sentinel: populated entries get overwritten below,
    // unpopulated ones keep -1.0 so the second pass can skip them via `>= 0`.
    double* hues = new double[_maxIter + 1ull];
    std::fill_n(hues, _maxIter + 1ull, -1.0);
    double hue_thr = 0.07;
    int total = 0;
    int item_count = 0;

    for (int i = 0; i < _maxIter; ++i) {
        total += _histogram[i];
    }

    // Entire view is inside the set: no histogram signal, fall back to a static palette.
    if (total == 0) {
        delete[] hues;
        CreateColorTables();
        return;
    }

    double hue = 0;
    for (int i = 0; i < _maxIter; ++i) {
        int item = _histogram[i];
        if (item == 0) {
            continue;
        }
        double d = (double)_histogram[i] / total;
        if (d > hue_thr)
            d = hue_thr;
        ++item_count;
        hues[i] = std::min(hue, 1.0);
        hue += d;
    }

    // Distribute the remaining hue budget (1.0 - hue) evenly across populated entries.
    // Each populated entry gets a cumulative correction of n*err, where n is its
    // ordinal among populated entries, so the final populated entry lands at ~1.0.
    double err = (1.0 - hue) / (item_count ? item_count : 1);
    int n = 0;
    for (int i = 0; i < _maxIter; ++i) {
        if (hues[i] >= 0) {
            ++n;
            hues[i] += err * n;
        }
    }

    // create HSV to QRgb table
    for (int i = 0; i < _maxIter; ++i) {
        double h = hues[i];
        if (h < 0)
            continue;
        h += offset;
        if (h >= 1)
            h -= 1;
        double section;
        double x = modf(h * 6, &section);
        int val = int(x * 255);
        switch (int(section) % 6) {
        case 0:
            _colorTable[i] = qRgb(255, 0, val);
            break;
        case 1:
            _colorTable[i] = qRgb(255 - val, 0, 255);
            break;
        case 2:
            _colorTable[i] = qRgb(0, val, 255);
            break;
        case 3:
            _colorTable[i] = qRgb(0, 255, 255 - val);
            break;
        case 4:
            _colorTable[i] = qRgb(val, 255, 0);
            break;
        case 5:
            _colorTable[i] = qRgb(255, 255 - val, 0);
            break;
        default:
            _colorTable[i] = qRgb(255, 255, 255);
        }
    }

    delete[] hues;
    setColorTableValid();
}

/**
 * @brief Build an iteration frequency histogram from the iteration buffer.
 *
 * Counts how many pixels escaped at each iteration count. Uses per-thread
 * private histograms with OpenMP to avoid contention, then merges them
 * in a critical section.
 *
 * @param pIterations Per-pixel iteration count buffer.
 * @param width Image width in pixels.
 * @param height Image height in pixels.
 */
void QMandelbrotWidget::CreateHistogram(const float* pIterations, int64_t width, int64_t height)
{
    delete[] _histogram;
    _histogram = nullptr;

    _histogram = new int[_maxIter + 1ull];
    memset(_histogram, 0, sizeof(int) * (_maxIter + 1ull));

#pragma omp parallel if (_useOpenMP)
    {
        int* histogram_private = new int[_maxIter + 1];
        memset(histogram_private, 0, sizeof(int) * (_maxIter + 1));

#pragma omp for
        for (int l = 0; l < height; ++l) {
            const float* pIter = pIterations + width * l;
            for (int k = 0; k < width; ++k) {
                int iter = (int)floorf(*pIter);
                iter = std::max(iter, 1);
                if (iter < _maxIter)
                    ++histogram_private[iter];
                ++pIter;
            }
        }
#pragma omp critical
        {
            for (int i = 0; i < _maxIter; ++i)
                _histogram[i] += histogram_private[i];
        }
        delete[] histogram_private;
    }
}

/**
 * @brief Convert the iteration count buffer to an RGB QImage.
 *
 * Maps each pixel's floating-point iteration count to a color via the
 * color table. Uses linear interpolation between adjacent color entries
 * for smooth coloring when fractional iteration counts are present.
 * Pixels that reached _maxIter (inside the set) are colored black.
 *
 * @param img Output QImage (must be Format_RGB32).
 * @param pIterations Per-pixel iteration count buffer.
 * @param width Image width in pixels.
 * @param height Image height in pixels.
 */
void QMandelbrotWidget::CreateDibFromIterations(QImage& img, const float* pIterations, int64_t width, int64_t height)
{
    // QImage is ARGB32 premultiplied?
    Q_ASSERT(img.format() == QImage::Format_RGB32);

#pragma omp parallel for schedule(static) if (_useOpenMP)
    for (int l = 0; l < height; ++l) {
        QRgb* scanLine = reinterpret_cast<QRgb*>(img.scanLine(l));
        const float* pIter = pIterations + width * l;
        for (int k = 0; k < width; ++k) {
            float mu = *pIter++;
            if (mu >= _maxIter) {
                scanLine[k] = qRgb(0, 0, 0);  // black
                continue;
            }
            float mu_i, mu_f = modff(mu, &mu_i);
            uint32_t index = (uint32_t)mu_i;
            if (index == (uint32_t)(_maxIter - 1)) {
                scanLine[k] = _colorTable[index];
            } else {
                QRgb c1 = _colorTable[index];
                QRgb c2 = _colorTable[index + 1];
                uint32_t alpha = (uint32_t)(256.0 * mu_f);
                scanLine[k] = blendAlphaQRgb(c1, c2, alpha);
            }
        }
    }
}

/**
 * @brief Render the fractal using IEEE 754 double-precision arithmetic.
 *
 * Implements the escape-time algorithm: Z(n+1) = Z(n)^2 + C, where
 * C is the pixel coordinate (Mandelbrot) or a fixed constant (Julia).
 * Iteration continues until |Z|^2 > 4 or the iteration limit is reached.
 *
 * Smooth coloring is computed as: mu = iter + 1 - log(log(|Z|)) / log(2).
 * X coordinates are pre-computed into a lookup table to avoid redundant
 * per-row calculation. Scanlines are parallelized via OpenMP.
 *
 * @param pIterations Output buffer for per-pixel iteration counts.
 * @param w Image width in pixels.
 * @param h Image height in pixels.
 * @param x0 Left edge of the view in the complex plane.
 * @param dx Horizontal step per pixel.
 * @param y0 Top edge of the view in the complex plane.
 * @param dy Vertical step per pixel.
 */
void QMandelbrotWidget::CalcIterationsDouble(float* pIterations, int64_t w, int64_t h, double x0, double dx, double y0, double dy)
{
    const float radius_sq = 2.0F * 2.0F;
    const float sqrt_32 = sqrt(32.f);
    bool isJulia = (_setType == stJulia);
    const double cr = isJulia ? _juliaConstant.real() : 0.0;
    const double ci = isJulia ? _juliaConstant.imag() : 0.0;

    double* xTable = new double[w];
    for (int i = 0; i < w; ++i) {
        xTable[i] = x0 + (double)i * dx;
    }

#pragma omp parallel for schedule(dynamic) if (_useOpenMP)
    for (int l = 0; l < h; ++l) {
        double y = y0 + (dy * l);
        double usq = 0, vsq = 0, u = 0, v = 0;
        double xc = (isJulia) ? cr : 0;
        double yc = (isJulia) ? ci : y;
        double modulus = 0;
        float* pbuff = pIterations + w * l;

        for (int k = 0; k < w; ++k) {
            int iter = 0;
            double x = xTable[k];
            if (isJulia) {
                u = x;
                v = y;
                usq = u * u;
                vsq = v * v;
                modulus = usq + vsq;
            } else {
                u = 0;
                v = 0;
                xc = x;
                usq = 0;
                vsq = 0;
                modulus = 0;
            }

            // Periodicity detection: inside-set orbits settle into a fixed cycle
            // in IEEE 754 precision. Snapshot (u, v) at exponentially spaced
            // iters and bail out when the current iterate matches the snapshot.
            double uRef = 0, vRef = 0;
            int period = 32;
            int nextSave = period;

            /*
                Complex iterative equation Z is:
                Mandebrot: Z(0) = 0, C = (x,y)
                Julia:     Z(0) = (x,y), C = Constant

                Shared:
                             2
                Z(i) = Z(i-1) + C

                check uv vector amplitude is smaller than 2
            */
            while (iter < _maxIter && modulus < radius_sq) {
                ++iter;

                // real
                double tmp = usq - vsq + xc;
                // imaginary:
                v = 2.0 * (u * v) + yc;
                u = tmp;
                vsq = v * v;
                usq = u * u;
                modulus = vsq + usq;

                if (u == uRef && v == vRef) {
                    iter = (int)_maxIter;
                    break;
                }
                if (iter == nextSave) {
                    uRef = u;
                    vRef = v;
                    period *= 2;
                    nextSave = iter + period;
                }
            }
            if (_smoothLevel && iter < _maxIter) {
                // modulus is in the range [4,36), create a scale between the 2 values.
                float mu = (float)(iter + 1) - ((float)sqrt(modulus - radius_sq)) / sqrt_32;
                *pbuff++ = mu;
            } else {
                *pbuff++ = (float)std::max(iter, 1);
            }
        }
    }

    delete[] xTable;
}

/**
 * @brief Render the fractal using 128-bit fixed-point precision.
 *
 * Same escape-time algorithm as CalcIterationsDouble() but uses fp128_t
 * for all complex-plane arithmetic, enabling extreme zoom levels
 * up to 2^113. The imaginary update uses a left-shift optimization:
 * v = (u * v) << 1 instead of v = 2 * u * v.
 *
 * @param pIterations Output buffer for per-pixel iteration counts.
 * @param width Image width in pixels.
 * @param height Image height in pixels.
 * @param x0 Left edge of the view in the complex plane.
 * @param dx Horizontal step per pixel.
 * @param y0 Top edge of the view in the complex plane.
 * @param dy Vertical step per pixel.
 */
void QMandelbrotWidget::CalcIterationsFP128(float* pIterations, int64_t width, int64_t height, fp128_t x0, fp128_t dx, fp128_t y0, fp128_t dy)
{
    const fp128_t radius_sq = 2 * 2;
    const float sqrt_32 = sqrt(32.f);
    bool isJulia = (_setType == stJulia);
    const fp128_t cr = isJulia ? _juliaConstant.real() : 0.0;
    const fp128_t ci = isJulia ? _juliaConstant.imag() : 0.0;

    fp128_t* xTable = new fp128_t[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

#pragma omp parallel for schedule(dynamic) if (_useOpenMP)
    for (int l = 0; l < height; ++l) {
        fp128_t y = y0 + (dy * l);
        fp128_t usq, vsq, u, v, x, tmp, modulus;
        fp128_t xc = (isJulia) ? cr : fp128_t(0);
        fp128_t yc = (isJulia) ? ci : y;
        float* pbuff = pIterations + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            x = xTable[k];
            // Julia
            if (isJulia) {
                u = x;
                v = y;
                usq = u * u;
                vsq = v * v;
                modulus = usq + vsq;
            }
            // Mandelbrot
            else {
                u = 0u;
                v = 0u;
                usq = 0u;
                vsq = 0u;
                xc = x;
                modulus = 0u;
            }

            // Periodicity detection: inside-set orbits collapse to a fixed cycle
            // in fp128 precision. Snapshot (u, v) at exponentially spaced iters
            // and bail out as soon as the current iterate matches the snapshot.
            // Cuts interior-pixel cost dramatically at deep zoom where _maxIter
            // is large.
            fp128_t uRef = 0u, vRef = 0u;
            int period = 32;
            int nextSave = period;

            /*
                Complex iterative equation Z is:
                Mandebrot: Z(0) = 0, C = (x,y)
                Julia:     Z(0) = (x,y), C = Constant

                Shared:
                             2
                Z(i) = Z(i-1) + C

                check uv vector amplitude is smaller than 2
            */
            while (iter < _maxIter && modulus < radius_sq) {
                ++iter;

                // real:
                tmp = usq - vsq + xc;
                // imaginary:
                // v = 2.0 * (u * v) + y;
                v = ((u * v) << 1) + yc;
                u = tmp;
                usq = u * u;
                vsq = v * v;
                modulus = usq + vsq;

                if (u == uRef && v == vRef) {
                    iter = (int)_maxIter;
                    break;
                }
                if (iter == nextSave) {
                    uRef = u;
                    vRef = v;
                    period *= 2;
                    nextSave = iter + period;
                }
            }

            if (_smoothLevel && iter < _maxIter) {
                // modulus is in the range [4,36), create a scale between the 2 values.
                float mu = (float)(iter + 1) -  (float)sqrt((float(modulus - radius_sq))) / sqrt_32;
                *pbuff++ = mu;
            } else {
                *pbuff++ = (float)std::max(iter, 1);
            }
        }
    }

    delete[] xTable;
}

/**
 * @brief Render the fractal and paint it to the widget surface.
 *
 * Allocates or reallocates the iteration buffer on resize, selects the
 * appropriate precision renderer based on zoom level, converts iteration
 * counts to an RGB image, and emits renderDone with frame statistics.
 *
 * @param event Paint event (unused).
 */
void QMandelbrotWidget::paintEvent(QPaintEvent* event)
{
    Q_UNUSED(event);
    QElapsedTimer paintTimer;
    paintTimer.start();

    QPainter p(this);
    if (_imageCache.isNull() || _imageCache.size() != size()) {
        _imageCache = QImage(size(), QImage::Format_RGB32);
        _imageCache.fill(Qt::white);

        delete[] _iterations;
        _iterations = nullptr;
        _iterations = new float[width() * height()];
        setFractalDataValid(false);
    }

    if (_autoIterations) {
        int64_t iters = calcAutoIterationLimits();
        if (iters != _maxIter) {
            _maxIter = iters;
            setFractalDataValid(false);
            setColorTableValid(false);
        }
    }

    int64_t w = width();
    int64_t h = height();
    bool historgamValid = (_paletteType == palHistogram) && fractalDataValid();

    if (!fractalDataValid()) {
        SetAspectRatio();
        fp128_t dx = (_xmax - _xmin) * (1.0 / w);
        fp128_t dy = dx;

        if (_precision == Precision::Double || (_precision == Precision::Auto && _zoomLevel <= MAX_DOUBLE_ZOOM_LEVEL)) {
            CalcIterationsDouble(_iterations, w, h, (double)_xmin, (double)dx, (double)_ymin, (double)dy);
        } else {
            CalcIterationsFP128(_iterations, w, h, _xmin, dx, _ymin, dy);
        }

        setFractalDataValid();
    }

    if (!historgamValid && _paletteType == palHistogram) {
        CreateHistogram(_iterations, w, h);
        CreateColorTableFromHistogram(_hsvOffset);
    } else if (!colorTableValid()) {
        CreateColorTables();
    }

    CreateDibFromIterations(_imageCache, _iterations, w, h);
    p.drawImage(rect().topLeft(), _imageCache);

    // prepare and emit frame stats
    FrameStats stats;
    qint64 elapsed = paintTimer.elapsed();
    stats.render_time_ms = static_cast<uint32_t>(elapsed > UINT32_MAX ? UINT32_MAX : elapsed);
    stats.zoom = static_cast<float>(_zoomLevel);
    stats.size = size();
    stats.max_iterations = _maxIter;
    emit renderDone(stats);
}

void QMandelbrotWidget::resizeEvent(QResizeEvent* event)
{
    QWidget::resizeEvent(event);
    invalidate(false);
}

/**
 * @brief Handle mouse clicks for zooming and view reset.
 *
 * Left click zooms in (2x, 4x with Ctrl, 8x with Ctrl+Shift).
 * Right click zooms out (2x, 4x with Ctrl, 8x with Ctrl+Shift).
 * Middle click resets to the default view.
 *
 * @param event Mouse event with button and modifier information.
 */
void QMandelbrotWidget::mousePressEvent(QMouseEvent* event)
{
    double zoomMultiplier = 0;
    if (event->button() == Qt::LeftButton) {
        zoomMultiplier = _zoomIncrement;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? _zoomIncrement * 4 : _zoomIncrement * 2;
        }
    } else if (event->button() == Qt::RightButton) {
        zoomMultiplier = 0.5;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? 1.0 / (_zoomIncrement * 4) : 1.0 / (_zoomIncrement * 2);
        }
    } else if (event->button() == Qt::MiddleButton) {
        resetView();
        return;
    }

    if (zoomMultiplier > 0) {
        _zoomLevel *= zoomMultiplier;
        if (_zoomLevel < 1.0) {
            _zoomLevel = 1.0;
        } else if (_zoomLevel > pow(2.0, logMaxZoom)) {
            _zoomLevel = pow(2.0, logMaxZoom);
        }

        OnZoomChange(event->pos(), zoomMultiplier);
    }
}

/**
 * @brief Compute automatic iteration limits based on zoom level.
 *
 * Linearly interpolates between min_iterations (at zoom 1x) and
 * max_iterations (at zoom 2^113) using the formula:
 * iters = min + (log2(zoom) / 113) * (max - min).
 *
 * @return The computed iteration limit.
 */
int64_t QMandelbrotWidget::calcAutoIterationLimits()
{
    // make iterations a function of zoom level.
    // map min_iterations to zoom=1 or smaller, and max_iterations to 2^113
    double logZoom = std::max(log2(_zoomLevel), 0.0);

    int64_t iters = static_cast<int64_t>(min_iterations + (logZoom / logMaxZoom) * (max_iterations - min_iterations));
    return iters;
}

/**
 * @brief Export the current view as a PNG file.
 *
 * Opens a file dialog for the user to choose the save location, then renders
 * the fractal at the requested resolution and writes it to disk. The X bounds
 * of the current view are preserved; Y bounds are recomputed to match the
 * target aspect ratio about the current Y center so the render is consistent
 * with what is displayed. The active palette (including histogram mode) is
 * honored.
 *
 * @param width Image width in pixels.
 * @param height Image height in pixels.
 */
void QMandelbrotWidget::saveImage(int width, int height)
{
    if (width <= 0 || height <= 0)
        return;

    QString fn = QFileDialog::getSaveFileName(this, tr("Save Image"), QString(), tr("PNG Files (*.png)"));
    if (fn.isEmpty())
        return;

    // Recompute Y bounds for the target aspect ratio around the current Y center.
    fp128_t ratio = (double)height / (double)width;
    fp128_t ysize = (_xmax - _xmin) * (ratio >> 1);
    fp128_t yCenter = (_ymax + _ymin) >> 1;
    fp128_t yMin = yCenter - ysize;

    fp128_t dx = (_xmax - _xmin) * (1.0 / width);
    fp128_t dy = dx;

    auto iterations = std::make_unique<float[]>((size_t)width * (size_t)height);

    if (_precision == Precision::Double || (_precision == Precision::Auto && _zoomLevel <= MAX_DOUBLE_ZOOM_LEVEL)) {
        CalcIterationsDouble(iterations.get(), width, height, (double)_xmin, (double)dx, (double)yMin, (double)dy);
    } else {
        CalcIterationsFP128(iterations.get(), width, height, _xmin, dx, yMin, dy);
    }

    if (_paletteType == palHistogram) {
        CreateHistogram(iterations.get(), width, height);
        CreateColorTableFromHistogram(_hsvOffset);
    }

    QImage img(width, height, QImage::Format_RGB32);
    CreateDibFromIterations(img, iterations.get(), width, height);
    img.save(fn, "PNG");

    // Histogram palette and color table were rebuilt against the export-resolution
    // iterations; trigger a recompute so the on-screen view is refreshed.
    if (_paletteType == palHistogram) {
        invalidate();
    }
}

void QMandelbrotWidget::setJuliaConstant(const std::complex<double>& c)
{
    SetDefaultValues();
    _juliaConstant = c;
    invalidate(false);
}

void QMandelbrotWidget::setSetType(set_type_t type)
{
    if (type >= stCount || type == _setType)
        return;

    _setType = type;
    resetView();
}

void QMandelbrotWidget::resetView()
{
    SetDefaultValues();
    invalidate();
}

void QMandelbrotWidget::zoomIn()
{
    _zoomLevel *= _zoomIncrement;
    OnZoomChange(QPoint(width() / 2, height() / 2), _zoomIncrement);
}

void QMandelbrotWidget::zoomOut()
{
    _zoomLevel /= _zoomIncrement;
    OnZoomChange(QPoint(width() / 2, height() / 2), 1.0 / _zoomIncrement);
}

/**
 * @brief Enable or disable palette color cycling animation.
 *
 * When enabled, starts a 30ms timer that rotates palette colors each tick.
 * When disabled, resets the color table to its static state and stops the timer.
 *
 * @param animate True to start animation, false to stop.
 */
void QMandelbrotWidget::setAnimatePalette(bool animate)
{
    _animate = animate;
    if (_animate) {
        // start timer
        _timer.setInterval(30ms);
        _timer.start();
    } else {
        // stop timer
        _hsvOffset = 0;
        _timer.stop();
    }
}

void QMandelbrotWidget::setPrecision(Precision p)
{
    _precision = p;
    invalidate();
}
void QMandelbrotWidget::setPaletteType(palette_t palette)
{
    _paletteType = palette;
    invalidate();
}

void QMandelbrotWidget::setSmoothTransitions(bool enable)
{
    if (enable != _smoothLevel) {
        _smoothLevel = enable;
        invalidate();
    }
}

/**
 * @brief Advance the palette animation by one frame.
 *
 * For histogram mode, rotates the HSV hue offset by 1/30.
 * For other modes, cyclically shifts all color table entries by one position.
 */
void QMandelbrotWidget::animationTick()
{
    // roll the _colorTable values

    if (_paletteType == palHistogram) {
        _hsvOffset += 1.0f / 30;
        CreateColorTableFromHistogram(_hsvOffset);
    } else {
        // Rotate visible palette entries [1 .. _maxIter - 1]. Leave entry 0
        // (white, escape-at-zero) and entry _maxIter (inside-set black) alone
        // so the hardcoded black does not leak into visible palette slots.
        QRgb first = _colorTable[1];
        for (int i = 1; i < _maxIter - 1; ++i) {
            _colorTable[i] = _colorTable[i + 1];
        }
        _colorTable[_maxIter - 1] = first;
    }
    update();
}

/**
 * @brief Pan the view horizontally by a fraction of the viewport width.
 * @param amount Fraction to pan (positive = right, negative = left).
 */
void QMandelbrotWidget::panHorizontal(double amount)
{
    auto dx = (_xmax - _xmin) * amount;
    _xmin += dx;
    _xmax += dx;
    invalidate(false);
}

/**
 * @brief Pan the view vertically by a fraction of the viewport height.
 * @param amount Fraction to pan (positive = down, negative = up).
 */
void QMandelbrotWidget::panVertical(double amount)
{
    auto dy = (_ymax - _ymin) * amount;
    _ymin += dy;
    _ymax += dy;
    invalidate(false);
}

/**
 * @brief Set the maximum iteration count or enable auto-iteration mode.
 *
 * A value of 0 enables automatic iteration scaling based on zoom level.
 * Any positive value sets a fixed iteration limit and rebuilds the color table.
 *
 * @param maxIter Iteration limit (0 = auto).
 */
void QMandelbrotWidget::setMaximumIterations(int64_t maxIter)
{
    if (maxIter < 0 || (maxIter == _maxIter && !_autoIterations))
        return;

    // Auto iterations
    _autoIterations = 0 == maxIter;
    if (_autoIterations) {
        invalidate();
        return;
    } else {
        _maxIter = maxIter;
    }

    invalidate();
}
