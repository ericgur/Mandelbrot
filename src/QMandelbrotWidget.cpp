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
#include <algorithm>
#include <cmath>

using namespace std::chrono_literals;

constexpr int32_t MAX_DOUBLE_LOG_ZOOM_LEVEL = 44;  ///< Deepest zoom, as a power of 2, that Auto precision renders with doubles.

/// Largest distance from the view center, on either axis, that PixelCoordinate() returns.
constexpr double pixelOffsetLimit = 5.0;
static_assert(QMandelbrotWidget::maxCenterMagnitude + pixelOffsetLimit < (1 << fp128IntBits), "every pixel coordinate must fit fp128_t");
static_assert(pixelOffsetLimit - QMandelbrotWidget::maxCenterMagnitude > 2, "a clamped pixel must escape at once");

/**
 * @brief Get the complex-plane coordinate of a pixel column or row.
 *
 * Pixel @p index of @p count sits at center + halfWidth * (2 * index - count) / width, which
 * puts pixel 0 on the view's edge and gives rows the same step as columns, 2 * halfWidth / width.
 *
 * The coordinate is accurate to a few LSBs at every zoom level. The ratio (2 * index - count) /
 * width splits into a whole number, which multiplies the half width exactly, and remainder /
 * width with |remainder| < width. Both operands of that division are scaled by the power of 2
 * just above width, which keeps them below 1 and exact in fp128_t, so the quotient and the
 * multiply after it are each rounded once, and the split keeps the fp128 factor in range in a
 * window of any shape. The alternatives are worse in two ways:
 * - Stepping from the edge by a rounded per-pixel step adds that step's rounding once per pixel,
 *   up to `width` LSBs across the image, which past zoom 2^100 shifts it by pixels and near
 *   logMaxZoom collapses the step to zero.
 * - Forming the ratio in double is off by 2^-53 of the half width, far above the LSB until zoom
 *   2^70 or so. That is invisible, but it decides the count of chaotic pixels near the boundary:
 *   at zoom 2^20, of the pixels where it disagreed with stepping, 10% matched an exact
 *   computation; with this division, all of them did.
 *
 * An offset beyond pixelOffsetLimit is clamped to it, which takes a window more than twice as
 * tall as it is wide at zoom 1x. Such a pixel is more than 3 from the origin on that axis
 * wherever the center is, clamped or not, so it escapes at once either way and only its smooth
 * shade moves. The clamp keeps every coordinate within fp128_t's range.
 *
 * @param center Center of the view on this axis.
 * @param halfWidth Half the view width.
 * @param index Pixel column or row.
 * @param count Image width for a column, image height for a row.
 * @param width Image width.
 * @return The coordinate on this axis.
 */
[[nodiscard]] static fp128_t PixelCoordinate(const fp128_t& center, const fp128_t& halfWidth, int64_t index, int64_t count, int64_t width)
{
    const int64_t numerator = 2 * index - count;
    const double offset = static_cast<double>(halfWidth) * static_cast<double>(numerator) / static_cast<double>(width);
    if (std::fabs(offset) > pixelOffsetLimit) {
        return center + fp128_t(std::copysign(pixelOffsetLimit, offset));
    }

    const int64_t whole = numerator / width;
    const int64_t remainder = numerator - whole * width;
    const double scale = std::exp2(std::ceil(std::log2(static_cast<double>(width))));
    const fp128_t fraction = fp128_t(static_cast<double>(remainder) / scale) / fp128_t(static_cast<double>(width) / scale);
    return center + halfWidth * whole + halfWidth * fraction;
}

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
 * @brief Initialize the view to the default complex plane region.
 *
 * Centers the view on the origin at zoom 1x, where it spans x in [-2.5, 2.5].
 */
void QMandelbrotWidget::SetDefaultValues()
{
    _logZoomLevel = 0;
    SetViewCenter(fp128_t(0), fp128_t(0));
}

void QMandelbrotWidget::SetViewCenter(const fp128_t& centerX, const fp128_t& centerY)
{
    constexpr fp128_t limit = maxCenterMagnitude;
    _centerX = std::clamp(centerX, -limit, limit);
    _centerY = std::clamp(centerY, -limit, limit);
}

fp128_t QMandelbrotWidget::ViewHalfWidth() const
{
    static_assert(logMinZoom >= 0, "the half-width is a right shift by the zoom level");
    static_assert(logMaxZoom + 1 <= fp128_t::F, "2.5 / 2^logMaxZoom must be exact in fp128_t");

    // The default view spans x in [-2.5, 2.5], and every zoom step halves that span. 2.5 is
    // 5 / 2, so the shifted value's lowest bit is 2^-(_logZoomLevel + 1) and the shift is exact.
    fp128_t halfWidth = 2.5;
    halfWidth >>= _logZoomLevel;
    return halfWidth;
}

/**
 * @brief Zoom by a power of 2 and center the view on a specific screen coordinate.
 *
 * The complex-plane coordinate of the given pixel becomes the new view center. The new zoom
 * level is clamped to [logMinZoom, logMaxZoom]; when the clamp leaves it unchanged, the view is
 * left as it is.
 *
 * @param point Screen coordinates of the new view center.
 * @param logZoomDelta Zoom steps of 2x each (positive to zoom in, negative to zoom out).
 */
void QMandelbrotWidget::OnZoomChange(const QPoint& point, int32_t logZoomDelta)
{
    const int32_t logZoomLevel = std::clamp(_logZoomLevel + logZoomDelta, logMinZoom, logMaxZoom);
    if (logZoomLevel == _logZoomLevel) {
        return;
    }

    const fp128_t halfWidth = ViewHalfWidth();
    const fp128_t centerX = PixelCoordinate(_centerX, halfWidth, point.x(), width(), width());
    const fp128_t centerY = PixelCoordinate(_centerY, halfWidth, point.y(), height(), width());

    _logZoomLevel = logZoomLevel;
    SetViewCenter(centerX, centerY);
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
    if (_setType == stJulia) {
        CalcIterationsDoubleImpl<true>(pIterations, w, h, x0, dx, y0, dy);
    } else {
        CalcIterationsDoubleImpl<false>(pIterations, w, h, x0, dx, y0, dy);
    }
}

template<bool IsJulia>
void QMandelbrotWidget::CalcIterationsDoubleImpl(float* pIterations, int64_t w, int64_t h, double x0, double dx, double y0, double dy)
{
    const float radius_sq = 2.0F * 2.0F;
    const float sqrt_32 = sqrt(32.f);
    const double cr = IsJulia ? static_cast<double>(_juliaConstant.real) : 0.0;
    const double ci = IsJulia ? static_cast<double>(_juliaConstant.imag) : 0.0;

    double* xTable = new double[w];
    for (int i = 0; i < w; ++i) {
        xTable[i] = x0 + (double)i * dx;
    }

#pragma omp parallel for schedule(dynamic) if (_useOpenMP)
    for (int l = 0; l < h; ++l) {
        const double y = y0 + (dy * l);
        const double yc = IsJulia ? ci : y;
        float* pbuff = pIterations + w * l;

        for (int k = 0; k < w; ++k) {
            int iter = 0;
            const double x = xTable[k];
            double u, v, usq, vsq, modulus, xc;

            if constexpr (IsJulia) {
                u = x;
                v = y;
                xc = cr;
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
                const double tmp = usq - vsq + xc;
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
 * up to 2^logMaxZoom. The imaginary update uses a left-shift optimization:
 * v = (u * v) << 1 instead of v = 2 * u * v.
 *
 * @param pIterations Output buffer for per-pixel iteration counts.
 * @param width Image width in pixels.
 * @param height Image height in pixels.
 * @param centerX Real part of the view center.
 * @param centerY Imaginary part of the view center.
 * @param halfWidth Half the view width in the complex plane.
 */
void QMandelbrotWidget::CalcIterationsFP128(float* pIterations, int64_t width, int64_t height, const fp128_t& centerX, const fp128_t& centerY,
                                            const fp128_t& halfWidth)
{
    if (_setType == stJulia) {
        CalcIterationsFP128Impl<true>(pIterations, width, height, centerX, centerY, halfWidth);
    } else {
        CalcIterationsFP128Impl<false>(pIterations, width, height, centerX, centerY, halfWidth);
    }
}

/**
 * @brief Test whether an iterate is still within the escape radius.
 *
 * The test is |Z|^2 < 4, but the iterate that escapes can have a squared modulus well beyond
 * fp128_t's range of [-2^fp128IntBits, 2^fp128IntBits): |Z| < 2 before a step and Julia constant
 * parts in [-2, 2] bound it by (4 + 2.83)^2 = 46.6, and a pixel PixelCoordinate() puts 7 from
 * the origin on both axes starts at 98. The multiply wraps silently, so that modulus can read
 * below 4 and keep a pixel iterating on garbage. u and v are therefore tested too: an iterate
 * with either part outside [-2, 2) has escaped whatever its modulus, and while both parts are
 * inside it the modulus is at most 8 - exact, or with 3 integer bits wrapped to -8, which reads
 * as escaped like the 8 it stands for. Where nothing overflowed this is the plain |Z|^2 < 4
 * test bit for bit, for every fp128IntBits from 3 up.
 *
 * All three tests read only the high QWORDs, which carry the sign and the integer bits; 2 and 4
 * have all-zero low QWORDs. Offsetting a part by 2 maps [-2, 2) onto [0, 4) as an unsigned
 * value, a modulus that hasn't wrapped is never negative, and 4 is a power of 2, so the three
 * values are all below 4 exactly when their OR is. That is two adds, two ORs and one compare,
 * and the fp128 render measured 1-3% faster with it than with the bare modulus compare it
 * replaced; four signed compares and their branches made it 7.5% slower.
 *
 * @param u Real part of the iterate.
 * @param v Imaginary part of the iterate.
 * @param modulus u^2 + v^2, as fp128 computed it.
 * @return True while the iterate has not escaped.
 */
[[nodiscard]] static FP128_FORCE_INLINE bool Bounded(const fp128_t& u, const fp128_t& v, const fp128_t& modulus)
{
    // 2 and 4 as the high QWORD of an fp128_t
    constexpr uint64_t two = 2ull << (fp128_t::F - 64);
    constexpr uint64_t four = 4ull << (fp128_t::F - 64);

    uint64_t low = 0, uHigh = 0, vHigh = 0, modulusHigh = 0;
    u.get_components(low, uHigh);
    v.get_components(low, vHigh);
    modulus.get_components(low, modulusHigh);

    return ((uHigh + two) | (vHigh + two) | modulusHigh) < four;
}

/**
 * @brief Get the smooth iteration count of an escaped iterate.
 *
 * mu = iter + 1 - sqrt(|Z|^2 - 4) / sqrt(32), with |Z|^2 in [4, 36) for Mandelbrot. The modulus
 * is formed in double from the parts, since the fp128 one can have wrapped (see Bounded()).
 *
 * @param iter Iteration at which the orbit escaped.
 * @param u Real part of the escaped iterate.
 * @param v Imaginary part of the escaped iterate.
 * @return The fractional iteration count.
 */
[[nodiscard]] static float SmoothIterations(int64_t iter, const fp128_t& u, const fp128_t& v)
{
    const double ud = static_cast<double>(u);
    const double vd = static_cast<double>(v);
    // rounding can take a modulus of exactly 4 a hair below it
    const double excess = std::max(ud * ud + vd * vd - 4.0, 0.0);
    return static_cast<float>(iter + 1) - std::sqrt(static_cast<float>(excess)) / std::sqrt(32.0f);
}

/**
 * @struct FP128Orbit
 * @brief Escape-time state of one pixel for the 128-bit fixed-point renderer.
 *
 * usq and vsq never leave an iteration: they are consumed only as their difference (the real
 * part of Z squared) and their sum (the escape test). The state therefore carries those two
 * instead, which is four QWORDs across the loop back edge rather than six. Both operations are
 * exact in fixed point, so the results are bit identical to computing them from usq and vsq.
 *
 * The members carry no default initializers on purpose: MakeOrbit() assigns every one of them,
 * and this type is constructed once per pixel on the hottest path in the renderer.
 */
struct FP128Orbit {
    fp128_t u, v;        ///< Current iterate Z = (u, v).
    fp128_t diff;        ///< usq - vsq, the real part of Z squared.
    fp128_t modulus;     ///< usq + vsq, compared against the escape radius.
    fp128_t uRef, vRef;  ///< Periodicity snapshot of an earlier iterate.
    int64_t iter;        ///< Iterations performed so far.
    int64_t period;      ///< Spacing to the next snapshot.
    int64_t nextSave;    ///< Iteration at which the next snapshot is taken.
    bool periodic;       ///< True once the orbit was found to repeat.
};

/// @brief Iterations run between periodicity checks. See RunOrbit() for why this is free to choose.
static constexpr int64_t fp128PeriodChunk = 16;

/**
 * @brief Pixels iterated in lockstep by the fixed-point renderer.
 *
 * The escape-time recurrence is a serial dependency chain, so a single pixel leaves most of the
 * machine idle no matter how few registers it needs. Stepping two neighbouring pixels together
 * gives the out-of-order engine two independent chains to overlap, which is worth ~19% to Clang
 * on a wide core. MSVC cannot schedule the doubled state and comes out 3-5% behind its own
 * single-pixel loop, so it keeps that one.
 */
static constexpr int fp128OrbitLanes =
#if defined(__clang__)
    2;
#else
    1;
#endif

/**
 * @brief Initialize an orbit for one pixel.
 * @tparam IsJulia True for Julia, false for Mandelbrot.
 * @param x Pixel real coordinate.
 * @param y Pixel imaginary coordinate.
 * @return The orbit at iteration zero.
 */
template<bool IsJulia>
[[nodiscard]] static FP128_FORCE_INLINE FP128Orbit MakeOrbit(const fp128_t& x, const fp128_t& y)
{
    FP128Orbit o;
    if constexpr (IsJulia) {
        // Julia starts at Z(0) = (x, y), so the first iteration needs its square already.
        o.u = x;
        o.v = y;
        const fp128_t usq = sqr(o.u);
        const fp128_t vsq = sqr(o.v);
        o.diff = usq - vsq;
        o.modulus = usq + vsq;
    } else {
        // Mandelbrot starts at Z(0) = 0, whose square is zero.
        o.u = 0u;
        o.v = 0u;
        o.diff = 0u;
        o.modulus = 0u;
    }

    o.uRef = 0u;
    o.vRef = 0u;
    o.iter = 0;
    o.period = 32;
    o.nextSave = o.period;
    o.periodic = false;

    return o;
}

/**
 * @brief Advance an orbit by one escape-time iteration.
 *
 * Z(i) = Z(i-1)^2 + C, with the imaginary part using a left shift for the doubling:
 * v = (u * v) << 1 instead of v = 2 * u * v.
 *
 * @param o Orbit to advance.
 * @param xc Real part of C.
 * @param yc Imaginary part of C.
 */
static FP128_FORCE_INLINE void StepOrbit(FP128Orbit& o, const fp128_t& xc, const fp128_t& yc)
{
    const fp128_t tmp = o.diff + xc;
    o.v = ((o.u * o.v) << 1) + yc;
    o.u = tmp;
    const fp128_t usq = sqr(o.u);
    const fp128_t vsq = sqr(o.v);
    o.diff = usq - vsq;
    o.modulus = usq + vsq;
}

/**
 * @brief Consult the periodicity snapshot and refresh it when due.
 *
 * Inside-set orbits collapse to an exact cycle in fp128 precision. Snapshots are taken at
 * exponentially spaced iterations, and an iterate equal to the snapshot means the orbit repeats
 * forever from there.
 *
 * @param o Orbit to test.
 * @return True when the orbit was found to repeat.
 */
static FP128_FORCE_INLINE bool CheckPeriodicity(FP128Orbit& o)
{
    if (o.u == o.uRef && o.v == o.vRef) {
        o.periodic = true;
        return true;
    }

    if (o.iter >= o.nextSave) {
        o.uRef = o.u;
        o.vRef = o.v;
        o.period *= 2;
        o.nextSave = o.iter + o.period;
    }

    return false;
}

/**
 * @brief Run one orbit to escape, to the iteration limit or to a detected cycle.
 *
 * The periodicity test is applied every fp128PeriodChunk iterations rather than every one, which
 * keeps the snapshot out of the innermost loop where it would occupy four of the sixteen general
 * registers. Detection is delayed by at most one chunk and that costs nothing: finding a cycle
 * only ever converts "this pixel will reach the iteration limit" into "this pixel reached it",
 * so a later detection produces the same value, just after a few more iterations.
 *
 * @param o Orbit to run.
 * @param xc Real part of C.
 * @param yc Imaginary part of C.
 * @param maxIter Iteration limit.
 */
static FP128_FORCE_INLINE void RunOrbit(FP128Orbit& o, const fp128_t& xc, const fp128_t& yc, int64_t maxIter)
{
    while (!o.periodic && o.iter < maxIter && Bounded(o.u, o.v, o.modulus)) {
        const int64_t stop = std::min(o.iter + fp128PeriodChunk, maxIter);

        while (o.iter < stop && Bounded(o.u, o.v, o.modulus)) {
            ++o.iter;
            StepOrbit(o, xc, yc);
        }

        if (CheckPeriodicity(o)) {
            break;
        }
    }
}

template<bool IsJulia>
void QMandelbrotWidget::CalcIterationsFP128Impl(float* pIterations, int64_t width, int64_t height, const fp128_t& centerX, const fp128_t& centerY,
                                                const fp128_t& halfWidth)
{
    const fp128_t cr = IsJulia ? _juliaConstant.real : fp128_t {};
    const fp128_t ci = IsJulia ? _juliaConstant.imag : fp128_t {};

    fp128_t* xTable = new fp128_t[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = PixelCoordinate(centerX, halfWidth, i, width, width);
    }

#pragma omp parallel for schedule(dynamic) if (_useOpenMP)
    for (int l = 0; l < height; ++l) {
        const fp128_t y = PixelCoordinate(centerY, halfWidth, l, height, width);
        const fp128_t yc = IsJulia ? ci : y;
        float* pbuff = pIterations + width * l;

        /*
            Complex iterative equation Z is:
            Mandebrot: Z(0) = 0, C = (x,y)
            Julia:     Z(0) = (x,y), C = Constant

            Shared:
                         2
            Z(i) = Z(i-1) + C

            check uv vector amplitude is smaller than 2
        */

        // Write one finished orbit out as a smooth or integer iteration count.
        // (not named "emit": Qt defines that as a macro for the signal keyword.)
        const auto writeResult = [&](const FP128Orbit& o) {
            const int64_t iter = o.periodic ? _maxIter : o.iter;
            if (_smoothLevel && iter < _maxIter) {
                *pbuff++ = SmoothIterations(iter, o.u, o.v);
            } else {
                *pbuff++ = (float)std::max<int64_t>(iter, 1);
            }
        };

        int k = 0;

        if constexpr (fp128OrbitLanes == 2) {
            // Two pixels in lockstep while both are still running; whichever outlives the other
            // finishes in the single-pixel loop below.
            for (; k + 1 < width; k += 2) {
                const fp128_t xcA = IsJulia ? cr : xTable[k];
                const fp128_t xcB = IsJulia ? cr : xTable[k + 1];
                FP128Orbit a = MakeOrbit<IsJulia>(xTable[k], y);
                FP128Orbit b = MakeOrbit<IsJulia>(xTable[k + 1], y);

                while (a.iter < _maxIter && Bounded(a.u, a.v, a.modulus) && b.iter < _maxIter && Bounded(b.u, b.v, b.modulus)) {
                    const int64_t stop = std::min(a.iter + fp128PeriodChunk, _maxIter);

                    // The two lanes advance together, so one iteration counter serves both.
                    while (a.iter < stop && Bounded(a.u, a.v, a.modulus) && Bounded(b.u, b.v, b.modulus)) {
                        ++a.iter;
                        ++b.iter;
                        StepOrbit(a, xcA, yc);
                        StepOrbit(b, xcB, yc);
                    }

                    if (CheckPeriodicity(a) || CheckPeriodicity(b)) {
                        break;
                    }
                }

                RunOrbit(a, xcA, yc, _maxIter);
                RunOrbit(b, xcB, yc, _maxIter);
                writeResult(a);
                writeResult(b);
            }
        }

        for (; k < width; ++k) {
            const fp128_t xc = IsJulia ? cr : xTable[k];
            FP128Orbit o = MakeOrbit<IsJulia>(xTable[k], y);

            RunOrbit(o, xc, yc, _maxIter);
            writeResult(o);
        }
    }

    delete[] xTable;
}

/**
 * @brief Compute a single Mandelbrot pixel at full fp128 precision.
 *
 * Used by the perturbation path as a fallback for glitched pixels and pixels
 * whose orbit outlives the reference. Mirrors the inner loop of
 * CalcIterationsFP128Impl<false>() but for one pixel.
 */
float QMandelbrotWidget::CalcSinglePixelFP128(fp128_t cx, fp128_t cy)
{
    fp128_t u = 0u, v = 0u, usq = 0u, vsq = 0u, modulus = 0u, tmp;
    fp128_t uRef = 0u, vRef = 0u;
    int period = 32;
    int nextSave = period;
    int iter = 0;

    while (iter < _maxIter && Bounded(u, v, modulus)) {
        ++iter;
        tmp = usq - vsq + cx;
        v = ((u * v) << 1) + cy;
        u = tmp;
        usq = sqr(u);
        vsq = sqr(v);
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
        return SmoothIterations(iter, u, v);
    }
    return (float)std::max(iter, 1);
}

/**
 * @brief Render the Mandelbrot set using perturbation theory.
 *
 * Computes one fp128 reference orbit at the view center, then iterates each
 * pixel as a small @c double delta against that orbit using:
 * @code
 *     Δ_{n+1} = 2·Z_n·Δ_n + Δ_n² + δc
 *     A_n     = Z_n + Δ_n   (actual orbit)
 * @endcode
 * Each iteration is 8 double multiplies — vs. 3 fp128 multiplies in the
 * standard path — so the speedup grows with the cost ratio of fp128:double
 * arithmetic (roughly 30–100× per multiply at this width).
 *
 * Pixels that satisfy Pauldelbrot's glitch criterion (|A|² ≪ |Z|²) or that
 * iterate past the reference's escape point fall back to per-pixel fp128.
 * Glitched-pixel fallback is correct but slow; for production use you'd
 * rebase the reference rather than rerun fp128, but this version is enough
 * for performance comparison.
 *
 * Julia is not handled by this path and falls through to CalcIterationsFP128Impl<true>().
 */
void QMandelbrotWidget::CalcIterationsPerturbation(float* pIterations, int64_t w, int64_t h, const fp128_t& centerX, const fp128_t& centerY,
                                                   const fp128_t& halfWidth)
{
    if (_setType == stJulia) {
        CalcIterationsFP128Impl<true>(pIterations, w, h, centerX, centerY, halfWidth);
        return;
    }

    const double radius_sq = 4.0;
    const float sqrt_32 = sqrt(32.f);
    constexpr double glitch_eps = 1e-6;

    // 1. Reference orbit at the view center, computed in fp128. Save the
    //    real/imag parts and modulus as doubles for the inner perturbation loop.
    const int64_t refX = w / 2;
    const int64_t refY = h / 2;
    const fp128_t cxRef = PixelCoordinate(centerX, halfWidth, refX, w, w);
    const fp128_t cyRef = PixelCoordinate(centerY, halfWidth, refY, h, w);

    auto refZx = std::make_unique<double[]>((size_t)_maxIter + 1);
    auto refZy = std::make_unique<double[]>((size_t)_maxIter + 1);
    auto refZmod = std::make_unique<double[]>((size_t)_maxIter + 1);
    refZx[0] = refZy[0] = refZmod[0] = 0.0;

    int64_t refLen = _maxIter;
    {
        fp128_t zx = 0u, zy = 0u, zxsq = 0u, zysq = 0u, zmod = 0u, ztmp;
        for (int64_t n = 1; n <= _maxIter; ++n) {
            ztmp = zxsq - zysq + cxRef;
            zy = ((zx * zy) << 1) + cyRef;
            zx = ztmp;
            zxsq = zx * zx;
            zysq = zy * zy;
            zmod = zxsq + zysq;

            const double zxd = (double)zx;
            const double zyd = (double)zy;
            refZx[n] = zxd;
            refZy[n] = zyd;
            refZmod[n] = zxd * zxd + zyd * zyd;

            if (!Bounded(zx, zy, zmod)) {
                refLen = n;
                break;
            }
        }
    }
    const bool refEscaped = (refLen < _maxIter);

    // 2. Each pixel's offset from the reference, halfWidth * 2 * (k - refX) / w, formed in double.
    //    halfWidth is a power of 2 times 5/2 and exact in double, so the offset is good to 53 bits
    //    at any zoom; it stays clear of fp128_t's narrow integer range, and differs from the
    //    difference of the fp128 coordinates only by their rounding.
    const double halfWidthD = static_cast<double>(halfWidth);
    auto xTable = std::make_unique<fp128_t[]>((size_t)w);
    auto dcxTable = std::make_unique<double[]>((size_t)w);
    for (int64_t k = 0; k < w; ++k) {
        xTable[k] = PixelCoordinate(centerX, halfWidth, k, w, w);
        dcxTable[k] = halfWidthD * (2.0 * static_cast<double>(k - refX) / static_cast<double>(w));
    }

#pragma omp parallel for schedule(dynamic) if (_useOpenMP)
    for (int l = 0; l < h; ++l) {
        const fp128_t y = PixelCoordinate(centerY, halfWidth, l, h, w);
        const double dcy = halfWidthD * (2.0 * static_cast<double>(l - refY) / static_cast<double>(w));
        float* pbuff = pIterations + w * l;

        for (int k = 0; k < w; ++k) {
            const double dcx = dcxTable[k];

            double Dx = 0, Dy = 0;
            int iter = 0;
            bool escaped = false;
            bool glitched = false;
            double mod = 0;

            while (iter < refLen) {
                const double zx_n = refZx[iter];
                const double zy_n = refZy[iter];

                // Δ_{n+1} = 2·Z_n·Δ_n + Δ_n² + δc
                const double newDx = 2.0 * (zx_n * Dx - zy_n * Dy) + (Dx * Dx - Dy * Dy) + dcx;
                const double newDy = 2.0 * (zx_n * Dy + zy_n * Dx) + 2.0 * Dx * Dy + dcy;
                Dx = newDx;
                Dy = newDy;
                ++iter;

                const double ax = refZx[iter] + Dx;
                const double ay = refZy[iter] + Dy;
                mod = ax * ax + ay * ay;

                if (mod >= radius_sq) {
                    escaped = true;
                    break;
                }
                // Pauldelbrot glitch check: |A|² ≪ |Z|² means catastrophic
                // cancellation when forming A = Z + Δ. Bail out and rerun
                // this pixel at full fp128 precision.
                if (mod < glitch_eps * refZmod[iter]) {
                    glitched = true;
                    break;
                }
            }

            // Reference exhausted before this pixel escaped: no high-precision
            // data left to perturb against. Fall back to fp128.
            const bool refExhausted = !escaped && !glitched && refEscaped && iter == refLen;

            float result;
            if (glitched || refExhausted) {
                result = CalcSinglePixelFP128(xTable[k], y);
            } else if (escaped) {
                if (_smoothLevel) {
                    result = (float)(iter + 1) - (float)sqrt(mod - radius_sq) / sqrt_32;
                } else {
                    result = (float)std::max(iter, 1);
                }
            } else {
                // Reached _maxIter without escape — inside the set.
                result = (float)_maxIter;
            }

            *pbuff++ = result;
        }
    }
}

void QMandelbrotWidget::CalcIterations(float* pIterations, int64_t width, int64_t height)
{
    const fp128_t halfWidth = ViewHalfWidth();

    if (_precision == Precision::Double || (_precision == Precision::Auto && _logZoomLevel <= MAX_DOUBLE_LOG_ZOOM_LEVEL)) {
        // the view's top left corner and pixel step, as PixelCoordinate() places them
        const double dx = 2.0 * static_cast<double>(halfWidth) / static_cast<double>(width);
        const double x0 = static_cast<double>(_centerX) - static_cast<double>(halfWidth);
        const double y0 = static_cast<double>(_centerY) - static_cast<double>(halfWidth) * static_cast<double>(height) / static_cast<double>(width);
        CalcIterationsDouble(pIterations, width, height, x0, dx, y0, dx);
    } else if (_precision == Precision::FixedPoint128) {
        CalcIterationsFP128(pIterations, width, height, _centerX, _centerY, halfWidth);
    } else {
        // Auto past 2^44 and explicit Perturbation both land here.
        // Julia falls through to fp128 inside CalcIterationsPerturbation.
        CalcIterationsPerturbation(pIterations, width, height, _centerX, _centerY, halfWidth);
    }
}

/**
 * @brief Render one frame into the image cache.
 *
 * Allocates or reallocates the iteration buffer on resize, selects the
 * appropriate precision renderer based on zoom level, and converts iteration
 * counts to an RGB image. Shared by paintEvent() and renderOffscreen().
 *
 * @return Wall-clock render time in nanoseconds.
 */
int64_t QMandelbrotWidget::RenderFrame()
{
    QElapsedTimer renderTimer;
    renderTimer.start();

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
        CalcIterations(_iterations, w, h);
        setFractalDataValid();
    }

    if (!historgamValid && _paletteType == palHistogram) {
        CreateHistogram(_iterations, w, h);
        CreateColorTableFromHistogram(_hsvOffset);
    } else if (!colorTableValid()) {
        CreateColorTables();
    }

    CreateDibFromIterations(_imageCache, _iterations, w, h);

    return renderTimer.nsecsElapsed();
}

/**
 * @brief Render the fractal and paint it to the widget surface.
 *
 * Renders through RenderFrame(), blits the result and emits renderDone with
 * frame statistics.
 *
 * @param event Paint event (unused).
 */
void QMandelbrotWidget::paintEvent(QPaintEvent* event)
{
    Q_UNUSED(event);
    QPainter p(this);

    const int64_t elapsedNs = RenderFrame();
    p.drawImage(rect().topLeft(), _imageCache);

    // prepare and emit frame stats
    FrameStats stats;
    const int64_t elapsedMs = elapsedNs / 1000000;
    stats.render_time_ms = static_cast<uint32_t>(elapsedMs > UINT32_MAX ? UINT32_MAX : elapsedMs);
    stats.log2Zoom = _logZoomLevel;
    stats.size = size();
    stats.max_iterations = _maxIter;
    emit renderDone(stats);
}

double QMandelbrotWidget::renderOffscreen()
{
    return static_cast<double>(RenderFrame()) / 1000000.0;
}

void QMandelbrotWidget::setView(const fp128_t& centerX, const fp128_t& centerY, int32_t log2Zoom)
{
    _logZoomLevel = std::clamp(log2Zoom, logMinZoom, logMaxZoom);
    SetViewCenter(centerX, centerY);
    invalidate();
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
    if (event->button() == Qt::MiddleButton) {
        resetView();
        return;
    }

    // one step zooms 2x; Ctrl makes it 4x and Ctrl+Shift 8x
    int32_t logZoomDelta = 1;
    if (event->modifiers() & Qt::ControlModifier) {
        logZoomDelta = (event->modifiers() & Qt::ShiftModifier) ? 3 : 2;
    }

    if (event->button() == Qt::LeftButton) {
        OnZoomChange(event->pos(), logZoomDelta);
    } else if (event->button() == Qt::RightButton) {
        OnZoomChange(event->pos(), -logZoomDelta);
    }
}

/**
 * @brief Compute automatic iteration limits based on zoom level.
 *
 * Linearly interpolates between min_iterations (at zoom 1x) and
 * max_iterations (at zoom 2^logMaxZoom) using the formula:
 * iters = min + (log2(zoom) / logMaxZoom) * (max - min).
 *
 * @return The computed iteration limit.
 */
int64_t QMandelbrotWidget::calcAutoIterationLimits()
{
    return autoIterationLimit(_logZoomLevel);
}

int64_t QMandelbrotWidget::autoIterationLimit(int32_t log2Zoom)
{
    // make iterations a function of zoom level.
    // map min_iterations to zoom=1 or smaller, and max_iterations to 2^logMaxZoom
    const double logZoom = std::max(log2Zoom, 0);

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

    auto iterations = std::make_unique<float[]>((size_t)width * (size_t)height);
    CalcIterations(iterations.get(), width, height);

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

Complex128 QMandelbrotWidget::defaultJuliaConstant()
{
    // parsed from text, since neither part has an exact binary form and a double would pin
    // them to the wrong digits beyond the 17th
    static const Complex128 constant {fp128_t("0.285"), fp128_t("0.01")};
    return constant;
}

void QMandelbrotWidget::setJuliaConstant(const Complex128& c)
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
    OnZoomChange(QPoint(width() / 2, height() / 2), 1);
}

void QMandelbrotWidget::zoomOut()
{
    OnZoomChange(QPoint(width() / 2, height() / 2), -1);
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
    // the view is 2 * halfWidth wide
    SetViewCenter(_centerX + (ViewHalfWidth() << 1) * amount, _centerY);
    invalidate(false);
}

/**
 * @brief Pan the view vertically by a fraction of the viewport height.
 * @param amount Fraction to pan (positive = down, negative = up).
 */
void QMandelbrotWidget::panVertical(double amount)
{
    // the view is 2 * halfWidth wide and height / width of that high
    const double viewWidths = amount * static_cast<double>(height()) / static_cast<double>(std::max(width(), 1));
    SetViewCenter(_centerX, _centerY + (ViewHalfWidth() << 1) * viewWidths);
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
