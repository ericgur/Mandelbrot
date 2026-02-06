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

using namespace std;

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
static inline QRgb blendAlphaQRgb(QRgb c1, QRgb c2, unsigned int alpha)
{
    unsigned int rb1 = (0x100 - alpha) * ((c1 & 0xFF00FF));
    unsigned int rb2 = alpha * ((c2 & 0xFF00FF));
    unsigned int g1 = (0x100 - alpha) * ((c1 & 0x00FF00));
    unsigned int g2 = alpha * ((c2 & 0x00FF00));
    unsigned int rb = (((rb1 + rb2) >> 8) & 0xFF00FF);
    unsigned int g = (((g1 + g2) >> 8) & 0x00FF00);
    return (QRgb)(rb | g);
}

QMandelbrotWidget::QMandelbrotWidget(QWidget* parent) : QWidget(parent), m_Timer(this)
{
    SetDefaultValues();
    CreateColorTables();
    connect(&m_Timer, &QChronoTimer::timeout, this, &QMandelbrotWidget::animationTick);
}

QMandelbrotWidget::~QMandelbrotWidget()
{
    delete[] m_Iterations;
    delete[] m_Histogram;
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
    m_Xmin = x_init[0];
    m_Xmax = x_init[1];
    m_Ymin = y_init[0];
    m_Ymax = y_init[1];
    m_ZoomLevel = zoom_init;
#else
    m_Xmax = 2.5;
    m_Xmin = -m_Xmax;
    m_Ymax = m_Ymin = 0;
    m_ZoomLevel = 1;
#endif
    SetAspectRatio();
}

/**
 * @brief Adjust Y bounds to preserve the correct aspect ratio.
 *
 * Recomputes m_Ymin and m_Ymax based on the current X bounds and
 * the widget's width/height ratio, keeping the Y center unchanged.
 */
void QMandelbrotWidget::SetAspectRatio()
{
    QSize s = size();

    // use m_Xmin, m_Xmax and m_rect to determine m_Ymin and m_Ymax
    // check if window created
    if (0 == s.height())
        return;

    fp128_t ratio = (double)s.height() / s.width();
    fp128_t ysize = (m_Xmax - m_Xmin) * (ratio >> 1);
    m_Ymin = ((m_Ymax + m_Ymin) >> 1) - ysize;
    m_Ymax = m_Ymin + (ysize << 1);
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
    fp128_t quarter = (m_Ymax - m_Ymin) * (1.0 / (zoomMultiplier * 2.0));
    fp128_t center = alpha * m_Ymax + (one - alpha) * m_Ymin;

    m_Ymin = center - quarter;
    m_Ymax = center + quarter;

    // fix x coords
    alpha = (double)(point.x()) / (s.width() - 1);
    quarter = (m_Xmax - m_Xmin) * (1.0 / (zoomMultiplier * 2.0));
    center = alpha * m_Xmax + (one - alpha) * m_Xmin;
    m_Xmin = center - quarter;
    m_Xmax = center + quarter;

    m_NeedToRecompute = true;
    update();
}

/**
 * @brief Build the color lookup table for the current palette type.
 *
 * Generates m_MaxIter + 1 color entries based on the active palette:
 * - @b Grey: Smooth greyscale from white to dark grey.
 * - @b Vivid: 6-segment HSV rainbow cycle.
 * - @b Gradient: Progressive RGB channel shifts (+3, +5, -3).
 *
 * Entry 0 is always white (escaped at iteration 0) and entry m_MaxIter
 * is always black (inside the set).
 */
void QMandelbrotWidget::CreateColorTables()
{
    m_ColorTable.clear();
    m_ColorTable.resize(m_MaxIter + 1);

    if (m_PaletteType == palGrey) {
        for (int64_t i = 1; i <= m_MaxIter; ++i) {
            int c = 255 - (int)(215.0f * (float)i / (float)m_MaxIter);
            m_ColorTable[i] = qRgb(c, c, c);
        }
    } else if (m_PaletteType == palVivid) {
        float step = (m_MaxIter <= 256) ? (13.0f / 256.f) : (13.0f / 256.f);
        for (int i = 0; i <= m_MaxIter; ++i) {
            float h = step * i;
            float x = (1.0f - fabs(fmodf(h, 2) - 1.0f));
            switch ((int)floorf(h) % 6) {
            case 0:  // 0-60
                m_ColorTable[i] = qRgb(255, int(x * 255), 0);
                break;
            case 1:  // 60-120
                m_ColorTable[i] = qRgb(int(x * 255), 255, 0);
                break;
            case 2:  // 120-180
                m_ColorTable[i] = qRgb(0, 255, int(x * 255));
                break;
            case 3:  // 180-240
                m_ColorTable[i] = qRgb(0, int(x * 255), 255);
                break;
            case 4:  // 240-300
                m_ColorTable[i] = qRgb(int(x * 255), 0, 255);
                break;
            case 5:  // 300-360
                m_ColorTable[i] = qRgb(255, 0, int(x * 255));
                break;
            }
        }
    } else if (m_PaletteType == palGradient) {
        unsigned int r = 0;
        unsigned int g = 20;
        unsigned int b = 255;

        for (int64_t i = 1; i <= m_MaxIter; ++i) {
            r = (r + 3) & 0xFF;
            g = (g + 5) & 0xFF;
            b = (b - 3) & 0xFF;
            m_ColorTable[i] = qRgb(b, g, r);
        }
    }

    m_ColorTable[0] = qRgb(255, 255, 255);
    m_ColorTable[m_MaxIter] = qRgb(0, 0, 0);
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
    double* hues = new double[m_MaxIter + 1ull];
    memset(hues, 0, sizeof(double) * (m_MaxIter + 1ull));
    double hue_thr = 0.07f;
    int total = 0;
    int item_count = 0;

    for (int i = 0; i < m_MaxIter; ++i) {
        total += m_Histogram[i];
    }

    double hue = 0;
    for (int i = 0; i < m_MaxIter; ++i) {
        int item = m_Histogram[i];
        if (item == 0) {
            hues[i] = -1.0f;
            continue;
        }
        double d = (double)m_Histogram[i] / total;
        if (d > hue_thr)
            d = hue_thr;
        ++item_count;
        hues[i] = min(hue, 1.0);
        hue += d;
    }

    double err = (1.0f - hue) / (item_count ? item_count : 1);
    for (int i = 0; i < m_MaxIter; ++i) {
        if (hues[i] >= 0)
            hues[i] += err * i;
    }

    // create HSV to QRgb table
    for (int i = 0; i < m_MaxIter; ++i) {
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
            m_ColorTable[i] = qRgb(255, 0, val);
            break;
        case 1:
            m_ColorTable[i] = qRgb(255 - val, 0, 255);
            break;
        case 2:
            m_ColorTable[i] = qRgb(0, val, 255);
            break;
        case 3:
            m_ColorTable[i] = qRgb(0, 255, 255 - val);
            break;
        case 4:
            m_ColorTable[i] = qRgb(val, 255, 0);
            break;
        case 5:
            m_ColorTable[i] = qRgb(255, 255 - val, 0);
            break;
        default:
            m_ColorTable[i] = qRgb(255, 255, 255);
        }
    }

    delete[] hues;
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
    if (m_Histogram != nullptr)
        delete[] m_Histogram;

    m_Histogram = new int[m_MaxIter + 1ull];
    memset(m_Histogram, 0, sizeof(int) * (m_MaxIter + 1ull));

#pragma omp parallel if (m_UseOpenMP)
    {
        int* histogram_private = new int[m_MaxIter + 1];
        if (histogram_private != nullptr) {
            memset(histogram_private, 0, sizeof(int) * (m_MaxIter + 1));
#pragma omp for
            for (int l = 0; l < height; ++l) {
                const float* pIter = pIterations + width * l;
                for (int k = 0; k < width; ++k) {
                    int iter = (int)floorf(*pIter);
                    iter = max(iter, 1);
                    if (iter < m_MaxIter)
                        ++histogram_private[iter];
                    ++pIter;
                }
            }
#pragma omp critical
            {
                for (int i = 0; i < m_MaxIter; ++i)
                    m_Histogram[i] += histogram_private[i];
            }
            delete[] histogram_private;
        }
    }
}

/**
 * @brief Convert the iteration count buffer to an RGB QImage.
 *
 * Maps each pixel's floating-point iteration count to a color via the
 * color table. Uses linear interpolation between adjacent color entries
 * for smooth coloring when fractional iteration counts are present.
 * Pixels that reached m_MaxIter (inside the set) are colored black.
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

    for (int l = 0; l < height; ++l) {
        QRgb* scanLine = reinterpret_cast<QRgb*>(img.scanLine(l));
        const float* pIter = pIterations + width * l;
        for (int k = 0; k < width; ++k) {
            float mu = *pIter++;
            if (mu >= m_MaxIter) {
                scanLine[k] = qRgb(0, 0, 0);  // black
                continue;
            }
            float mu_i, mu_f = modff(mu, &mu_i);
            unsigned int index = (unsigned int)mu_i;
            if (index < 0) {  // unsigned cannot be <0, keep original logic
                scanLine[k] = m_ColorTable[0];
                continue;
            }
            if (index == (unsigned int)(m_MaxIter - 1)) {
                scanLine[k] = m_ColorTable[index];
            } else {
                QRgb c1 = m_ColorTable[index];
                QRgb c2 = m_ColorTable[index + 1];
                unsigned int alpha = (unsigned int)(256.0 * mu_f);
                if (alpha > 255) {
                    scanLine[k] = c2;
                } else {
                    scanLine[k] = blendAlphaQRgb(c1, c2, alpha);
                }
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
void QMandelbrotWidget::DrawImageDouble(float* pIterations, int64_t w, int64_t h, double x0, double dx, double y0, double dy)
{
    const float radius_sq = 2.0F * 2.0F;
    const float LOG2 = logf(2.0F);
    bool isJulia = (m_SetType == stJulia);
    const double cr = isJulia ? m_JuliaConstant.real() : 0.0;
    const double ci = isJulia ? m_JuliaConstant.imag() : 0.0;

    double* xTable = new double[w];
    for (int i = 0; i < w; ++i) {
        xTable[i] = x0 + (double)i * dx;
    }

#pragma omp parallel for schedule(dynamic) if (m_UseOpenMP)
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

            /*
                Complex iterative equation Z is:
                Mandebrot: Z(0) = 0, C = (x,y)
                Julia:     Z(0) = (x,y), C = Constant

                Shared:
                             2
                Z(i) = Z(i-1) + C

                check uv vector amplitude is smaller than 2
            */
            while (modulus < radius_sq && ++iter < m_MaxIter) {
                // real
                double tmp = usq - vsq + xc;
                // imaginary:
                // v = 2.0 * (u * v) + y;
                v = u * v + u * v + yc;
                u = tmp;
                vsq = v * v;
                usq = u * u;
                modulus = vsq + usq;
            }
            if (m_SmoothLevel && iter < m_MaxIter && iter > 0) {
                float mu = (float)(iter + 1) - (logf(logf(sqrtf((float)modulus)))) / LOG2;
                *pbuff++ = max(mu, 1.0f);
            } else {
                *pbuff++ = (float)max(iter, 1);
            }
        }
    }

    delete[] xTable;
}

/**
 * @brief Render the fractal using 128-bit fixed-point precision.
 *
 * Same escape-time algorithm as DrawImageDouble() but uses fp128_t
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
void QMandelbrotWidget::DrawImageFixedPoint128(float* pIterations, int64_t width, int64_t height, fp128_t x0, fp128_t dx, fp128_t y0, fp128_t dy)
{
    const fp128_t radius_sq = 2 * 2;
    const float LOG2 = logf(2.0F);
    bool isJulia = (m_SetType == stJulia);
    const fp128_t cr = isJulia ? m_JuliaConstant.real() : 0.0;
    const fp128_t ci = isJulia ? m_JuliaConstant.imag() : 0.0;

    fp128_t* xTable = new fp128_t[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

#pragma omp parallel for schedule(dynamic) if (m_UseOpenMP)
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

            /*
                Complex iterative equation Z is:
                Mandebrot: Z(0) = 0, C = (x,y)
                Julia:     Z(0) = (x,y), C = Constant

                Shared:
                             2
                Z(i) = Z(i-1) + C

                check uv vector amplitude is smaller than 2
            */
            while (modulus < radius_sq && ++iter < m_MaxIter) {
                // real:
                tmp = usq - vsq + xc;
                // imaginary:
                // v = 2.0 * (u * v) + y;
                v = ((u * v) << 1) + yc;
                u = tmp;
                usq = u * u;
                vsq = v * v;
                modulus = usq + vsq;
            }

            if (m_SmoothLevel && iter < m_MaxIter && iter > 1) {
                *pbuff++ = (float)(iter + 1) - (logf(logf(sqrtf((float)modulus)))) / LOG2;
            } else {
                *pbuff++ = (float)max(iter, 1);
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
    QElapsedTimer _paintTimer;
    _paintTimer.start();

    QPainter p(this);
    if (m_ImageCache.isNull() || m_ImageCache.size() != size()) {
        m_ImageCache = QImage(size(), QImage::Format_RGB32);
        m_ImageCache.fill(Qt::white);

        delete[] m_Iterations;
        m_Iterations = new float[width() * height()];
        m_NeedToRecompute = true;
    }

    if (m_AutoIterations) {
        int64_t iters = calcAutoIterationLimits();
        if (iters != m_MaxIter) {
            m_MaxIter = iters;
            CreateColorTables();
            m_NeedToRecompute = true;
        }
    }

    int64_t w = width();
    int64_t h = height();

    if (m_NeedToRecompute) {
        SetAspectRatio();
        fp128_t dx = (m_Xmax - m_Xmin) * (1.0 / w);
        fp128_t dy = dx;

        if (m_Precision == Precision::Double || (m_Precision == Precision::Auto && m_ZoomLevel <= (1ull << 44))) {
            DrawImageDouble(m_Iterations, w, h, (double)m_Xmin, (double)dx, (double)m_Ymin, (double)dy);
        } else {
            DrawImageFixedPoint128(m_Iterations, w, h, m_Xmin, dx, m_Ymin, dy);
        }

        if (m_PaletteType == palHistogram) {
            CreateHistogram(m_Iterations, w, h);
            CreateColorTableFromHistogram(m_HsvOffset);
        }

        m_NeedToRecompute = false;
    }

    CreateDibFromIterations(m_ImageCache, m_Iterations, w, h);
    p.drawImage(rect().topLeft(), m_ImageCache);

    // prepare and emit frame stats
    FrameStats stats;
    qint64 elapsed = _paintTimer.elapsed();
    stats.render_time_ms = static_cast<uint32_t>(elapsed > UINT32_MAX ? UINT32_MAX : elapsed);
    stats.zoom = static_cast<float>(m_ZoomLevel);
    stats.size = size();
    stats.max_iterations = m_MaxIter;
    emit renderDone(stats);
}

void QMandelbrotWidget::resizeEvent(QResizeEvent* event)
{
    QWidget::resizeEvent(event);
    m_NeedToRecompute = true;
    update();
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
    if (event->button() == Qt::LeftButton) {
        double zoomMultiplier = 2.0;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? 8.0 : 4.0;
        }
        m_ZoomLevel *= zoomMultiplier;
        OnZoomChange(event->pos(), zoomMultiplier);
    } else if (event->button() == Qt::RightButton) {
        double zoomMultiplier = 0.5;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? 0.125 : 0.25;
        }
        m_ZoomLevel *= zoomMultiplier;
        OnZoomChange(event->pos(), zoomMultiplier);
    } else if (event->button() == Qt::MiddleButton) {
        resetView();
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
    // map 64 iterations to zoom=1 or smaller, and max_iterations to 2^113
    double logZoom = max(log2(m_ZoomLevel), 0);

    int64_t iters = static_cast<int64_t>(min_iterations + (logZoom / logMaxZoom) * (max_iterations - min_iterations));
    return iters;
}

/**
 * @brief Export the current view as a PNG file.
 *
 * Opens a file dialog for the user to choose the save location,
 * then renders and saves the image at the specified resolution.
 *
 * @param width Image width in pixels.
 * @param height Image height in pixels.
 */
void QMandelbrotWidget::saveImage(int width, int height)
{
    QImage img(width, height, QImage::Format_ARGB32);
    img.fill(Qt::white);
    QString fn = QFileDialog::getSaveFileName(this, tr("Save Image"), QString(), tr("PNG Files (*.png)"));
    if (fn.isEmpty())
        return;
    img.save(fn, "PNG");
}

void QMandelbrotWidget::setJuliaConstant(const std::complex<double>& c)
{
    SetDefaultValues();
    m_JuliaConstant = c;
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::setSetType(set_type_t type)
{
    if (type >= stCount || type == m_SetType)
        return;

    m_SetType = type;
    resetView();
}

void QMandelbrotWidget::resetView()
{
    SetDefaultValues();
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::zoomIn()
{
    m_ZoomLevel *= 2.0;
    OnZoomChange(QPoint(width() / 2, height() / 2), m_ZoomIncrement);
}

void QMandelbrotWidget::zoomOut()
{
    m_ZoomLevel /= 2.0;
    OnZoomChange(QPoint(width() / 2, height() / 2), 1.0 / m_ZoomIncrement);
}

/**
 * @brief Enable or disable palette color cycling animation.
 *
 * When enabled, starts a 70ms timer that rotates palette colors each tick.
 * When disabled, resets the color table to its static state and stops the timer.
 *
 * @param animate True to start animation, false to stop.
 */
void QMandelbrotWidget::setAnimatePalette(bool animate)
{
    m_Animate = animate;
    if (m_Animate) {
        // start timer
        m_Timer.setInterval(70ms);
        m_Timer.start();
    } else {
        if (m_PaletteType == palHistogram) {
            m_HsvOffset = 0;
            CreateColorTableFromHistogram(m_HsvOffset);
        } else {
            CreateColorTables();  // reset color tables
        }

        m_Timer.stop();
    }

    update();
}

void QMandelbrotWidget::setPrecision(Precision p)
{
    m_Precision = p;
    m_NeedToRecompute = true;
    update();
}

/**
 * @brief Advance the palette animation by one frame.
 *
 * For histogram mode, rotates the HSV hue offset by 1/30.
 * For other modes, cyclically shifts all color table entries by one position.
 */
void QMandelbrotWidget::animationTick()
{
    // roll the m_ColorTable values

    if (m_PaletteType == palHistogram) {
        m_HsvOffset += 1.0f / 30;
        CreateColorTableFromHistogram(m_HsvOffset);
    } else {
        QVector<QRgb> temp = m_ColorTable;
        temp[0] = m_ColorTable[0];
        for (int i = 1; i < m_MaxIter; ++i) {
            temp[i] = m_ColorTable[i + 1];
        }
        temp[m_MaxIter] = m_ColorTable[1];

        m_ColorTable = std::move(temp);
    }
    repaint();
}

/**
 * @brief Pan the view horizontally by a fraction of the viewport width.
 * @param amount Fraction to pan (positive = right, negative = left).
 */
void QMandelbrotWidget::panHorizontal(double amount)
{
    auto dx = (m_Xmax - m_Xmin) * amount;
    m_Xmin += dx;
    m_Xmax += dx;
    m_NeedToRecompute = true;
    update();
}

/**
 * @brief Pan the view vertically by a fraction of the viewport height.
 * @param amount Fraction to pan (positive = down, negative = up).
 */
void QMandelbrotWidget::panVertical(double amount)
{
    auto dy = (m_Ymax - m_Ymin) * amount;
    m_Ymin += dy;
    m_Ymax += dy;
    m_NeedToRecompute = true;
    update();
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
    if (maxIter < 0 || maxIter == m_MaxIter)
        return;

    // Auto iterations
    m_AutoIterations = 0 == maxIter;
    if (m_AutoIterations) {
        update();
        return;
    } else {
        m_MaxIter = maxIter;
    }

    CreateColorTables();
    m_NeedToRecompute = true;
    update();
}
