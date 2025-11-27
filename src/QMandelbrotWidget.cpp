#include "pch.h"
#include "QMandelbrotWidget.h"
#include <QPainter>
#include <QFileDialog>
#include <QMouseEvent>
#include <QElapsedTimer>
#include <cmath>

using namespace std;

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

QMandelbrotWidget::QMandelbrotWidget(QWidget* parent)
    : QWidget(parent),
    m_Timer(this)
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

void QMandelbrotWidget::SetAspectRatio()
{
    QSize s = size();

    // use m_Xmin, m_Xmax and m_rect to determine m_Ymin and m_Ymax
    // check if window created
    if (0 == s.height())
        return;

    fixed_8_120_t ratio = (double)s.height() / s.width();
    fixed_8_120_t ysize = (m_Xmax - m_Xmin) * (ratio >> 1);
    m_Ymin = ((m_Ymax + m_Ymin) >> 1) - ysize;
    m_Ymax = m_Ymin + (ysize << 1);
}

void QMandelbrotWidget::OnZoomChange(const QPoint& point, double zoomMultiplier)
{
    QSize s = size();
    static fixed_8_120_t one = 1;

    // fix y coords
    fixed_8_120_t alpha = (double)(point.y()) / (s.height() - 1);
    fixed_8_120_t quarter = (m_Ymax - m_Ymin) * (1.0 / (zoomMultiplier * 2.0));
    fixed_8_120_t center = alpha * m_Ymax + (one - alpha) * m_Ymin;

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

void QMandelbrotWidget::CreateColorTables()
{
    m_ColorTable.clear();
    m_ColorTable.resize(m_MaxIter + 1);

    if (m_PaletteType == palGrey) {
        for (int64_t i = 1; i <= m_MaxIter; ++i) {
            int c = 255 - (int)(215.0f * (float)i / (float)m_MaxIter);
            m_ColorTable[i] = qRgb(c, c, c);
        }
    }
    else if (m_PaletteType == palVivid) {
        float step = (m_MaxIter <= 256) ? (13.0f / 256.f) : (13.0f / 256.f);
        for (int i = 0; i <= m_MaxIter; ++i) {
            float h = step * i;
            float x = (1.0f - fabs(fmodf(h, 2) - 1.0f));
            switch ((int)floorf(h) % 6) {
            case 0: // 0-60 
                m_ColorTable[i] = qRgb(255, int(x * 255), 0);
                break;
            case 1: // 60-120
                m_ColorTable[i] = qRgb(int(x * 255), 255, 0);
                break;
            case 2: // 120-180 
                m_ColorTable[i] = qRgb(0, 255, int(x * 255));
                break;
            case 3: // 180-240
                m_ColorTable[i] = qRgb(0, int(x * 255), 255);
                break;
            case 4: // 240-300
                m_ColorTable[i] = qRgb(int(x * 255), 0, 255);
                break;
            case 5: // 300-360
                m_ColorTable[i] = qRgb(255, 0, int(x * 255));
                break;
            }
        }
    }
    else if (m_PaletteType == palGradient) {
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

void QMandelbrotWidget::CreateHistogram(const float* pIterations, int64_t width, int64_t height)
{
    if (m_Histogram != nullptr)
        delete[] m_Histogram;

    m_Histogram = new int[m_MaxIter + 1ull];
    memset(m_Histogram, 0, sizeof(int) * (m_MaxIter + 1ull));

#pragma omp parallel 
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
                scanLine[k] = qRgb(0, 0, 0); // black
                continue;
            }
            float mu_i, mu_f = modff(mu, &mu_i);
            unsigned int index = (unsigned int)mu_i;
            if (index < 0) { // unsigned cannot be <0, keep original logic
                scanLine[k] = m_ColorTable[0];
                continue;
            }
            if (index == (unsigned int)(m_MaxIter - 1)) {
                scanLine[k] = m_ColorTable[index];
            }
            else {
                QRgb c1 = m_ColorTable[index];
                QRgb c2 = m_ColorTable[index + 1];
                unsigned int alpha = (unsigned int)(256.0 * mu_f);
                if (alpha > 255) {
                    scanLine[k] = c2;
                }
                else {
                    scanLine[k] = blendAlphaQRgb(c1, c2, alpha);
                }
            }
        }
    }
}

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

    #pragma omp parallel for schedule(dynamic)
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
            }
            else {
                u = 0; 
                v = 0; 
                xc = x; 
                usq = 0; 
                vsq = 0; 
                modulus = 0; 
            }

            // complex iterative equation is:
            // Z(i) = Z(i-1) ^ 2 + C
            // Mandebrot: Z(0) = 0, C = (x,y)
            // Julia:     Z(0) = (x,y), C = Constant
            // 
            // check uv vector amplitude is smaller than 2

            while (modulus < radius_sq && ++iter < m_MaxIter) {
                // real
                double tmp = usq - vsq + xc;
                // imaginary
                //v = 2.0 * (u * v) + y;
                v = u * v + u * v + yc;
                u = tmp;
                vsq = v * v;
                usq = u * u;
                modulus = vsq + usq;
            }
            if (m_SmoothLevel && iter < m_MaxIter && iter > 0) {
                float mu = (float)(iter + 1) - (logf(logf(sqrtf((float)modulus)))) / LOG2;
                *pbuff++ = max(mu, 1.0f);
            }
            else {
                *pbuff++ = (float)max(iter, 1);
            }
        }
    }

    delete[] xTable;
}

void QMandelbrotWidget::DrawImageFixedPoint128(float* pIterations, int64_t width, int64_t height, const fixed_8_120_t& x0, const fixed_8_120_t& dx,
    const fixed_8_120_t& y0, const fixed_8_120_t& dy)
{
    const fixed_8_120_t radius_sq = 2 * 2;
    const float LOG2 = logf(2.0F);
    bool isJulia = (m_SetType == stJulia);
    const fixed_8_120_t cr = isJulia ? m_JuliaConstant.real() : 0.0;
    const fixed_8_120_t ci = isJulia ? m_JuliaConstant.imag() : 0.0;

    fixed_8_120_t* xTable = new fixed_8_120_t[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

    #pragma omp parallel for schedule(dynamic)
    for (int l = 0; l < height; ++l) {
        fixed_8_120_t y = y0 + (dy * l);
        fixed_8_120_t usq, vsq, u, v, x, tmp, modulus;
        fixed_8_120_t xc = (isJulia) ? cr : fixed_8_120_t(0);
        fixed_8_120_t yc = (isJulia) ? ci : y;

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

            // complex iterative equation is:
            // Z(i) = Z(i-1) ^ 2 + C
            // Mandebrot: Z(0) = 0, C = (x,y)
            // Julia:     Z(0) = (x,y), C = Constant
            // 
            // check uv vector amplitude is smaller than 2

            while (modulus < radius_sq && ++iter < m_MaxIter) {
                // real
                tmp = usq - vsq + xc;
                // imaginary
                //v = 2.0 * (u * v) + y;
                v = ((u * v) << 1) + yc;
                u = tmp;
                usq = u * u;
                vsq = v * v;
                modulus = usq + vsq;
            }

            if (m_SmoothLevel && iter < m_MaxIter && iter > 1) {
                *pbuff++ = (float)(iter + 1) - (logf(logf(sqrtf((float)modulus)))) / LOG2;
            }
            else {
                *pbuff++ = (float)max(iter, 1);
            }
        }
    }

    delete[] xTable;
}

void QMandelbrotWidget::paintEvent(QPaintEvent* event)
{
    Q_UNUSED(event);
    QElapsedTimer _paintTimer; _paintTimer.start();

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
        fixed_8_120_t dx = (m_Xmax - m_Xmin) * (1.0 / w);
        fixed_8_120_t dy = dx;

        if (m_Precision == Precision::Double || (m_Precision == Precision::Auto && m_ZoomLevel <= (1ull<<44))) {
            DrawImageDouble(m_Iterations, w, h, (double)m_Xmin, (double)dx, (double)m_Ymin, (double)dy);
        }
        else {
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

void QMandelbrotWidget::mousePressEvent(QMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        double zoomMultiplier = 2.0;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? 8.0 : 4.0;
        }
        m_ZoomLevel *= zoomMultiplier;
        OnZoomChange(event->pos(), zoomMultiplier);
    }
    else if (event->button() == Qt::RightButton) {
        double zoomMultiplier = 0.5;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? 0.125 : 0.25;
        }
        m_ZoomLevel *= zoomMultiplier;
        OnZoomChange(event->pos(), zoomMultiplier);
    }
    else if (event->button() == Qt::MiddleButton) {
        resetView();
    }
}

int64_t QMandelbrotWidget::calcAutoIterationLimits()
{
    // make iterations a function of zoom level. 
    // map 64 iterations to zoom=1 or smaller, and max_iterations to 2^113
    double logZoom = max(log2(m_ZoomLevel), 0);
    
    int64_t iters = static_cast<int64_t>(min_iterations + (logZoom / logMaxZoom) * (max_iterations - min_iterations));
    return iters;
}

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
    OnZoomChange(QPoint(width() / 2, height() / 2), m_ZoomIncrement);
}

void QMandelbrotWidget::zoomOut()
{
    OnZoomChange(QPoint(width() / 2, height() / 2), 1.0 / m_ZoomIncrement);
}

void QMandelbrotWidget::setAnimatePalette(bool animate)
{
    m_Animate = animate;
    if (m_Animate) {
        // start timer
        m_Timer.setInterval(70ms);
        m_Timer.start();
    }
    else {
        if (m_PaletteType == palHistogram) {
            m_HsvOffset = 0;
            CreateColorTableFromHistogram(m_HsvOffset);
        }
        else {
            CreateColorTables(); // reset color tables
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

void QMandelbrotWidget::animationTick()
{
    // roll the m_ColorTable values

    if (m_PaletteType == palHistogram) {
        m_HsvOffset += 1.0f / 30;
        CreateColorTableFromHistogram(m_HsvOffset);
    }
    else {
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

void QMandelbrotWidget::panHorizontal(double amount)
{
    auto dx = (m_Xmax - m_Xmin) * amount;
    m_Xmin += dx;
    m_Xmax += dx;
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::panVertical(double amount)
{
    auto dy = (m_Ymax - m_Ymin) * amount;
    m_Ymin += dy;
    m_Ymax += dy;
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::setMaximumIterations(int64_t maxIter)
{
    if (maxIter < 0 || maxIter == m_MaxIter)
        return;

    // Auto iterations
    m_AutoIterations = 0 == maxIter;
    if (m_AutoIterations) {
        update();
        return; 
    }
    else {
        m_MaxIter = maxIter;
    }

    CreateColorTables();
    m_NeedToRecompute = true;
    update();
}
