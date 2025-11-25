#include "pch.h"
#include "QMandelbrotWidget.h"
#include <QPainter>
#include <QFileDialog>
#include <QImageWriter>
#include <QMouseEvent>
#include <QElapsedTimer>
#include <cmath>
#include <algorithm>

using namespace std;

static inline constexpr QRgb BGR_to_QRgb(int b, int g, int r)
{
    return qRgb(r, g, b);
}

static inline QRgb blendAlphaQRgb(QRgb colora, QRgb colorb, unsigned int alpha)
{
    unsigned int rb1 = (0x100 - alpha) * ((colora & 0xFF00FF));
    unsigned int rb2 = alpha * ((colorb & 0xFF00FF));
    unsigned int g1 = (0x100 - alpha) * ((colora & 0x00FF00));
    unsigned int g2 = alpha * ((colorb & 0x00FF00));
    unsigned int rb = (((rb1 + rb2) >> 8) & 0xFF00FF);
    unsigned int g = (((g1 + g2) >> 8) & 0x00FF00);
    return (QRgb)(rb | g);
}

QMandelbrotWidget::QMandelbrotWidget(QWidget* parent)
    : QWidget(parent),
    m_Timer(this)
{
    //setAutoFillBackground(true);
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
    m_xmin = x_init[0];
    m_xmax = x_init[1];
    m_ymin = y_init[0];
    m_ymax = y_init[1];
    m_zoom = zoom_init;
#else
    m_xmax = 2.5;
    m_xmin = -m_xmax;
    m_ymax = m_ymin = 0;
    m_zoom = 1;
#endif
    SetAspectRatio();
}

void QMandelbrotWidget::SetAspectRatio()
{
    QRect r = rect();

    //use m_xmin, m_xmax and m_rect to determine m_ymin and m_ymax
    //check if window created
    if (0 == r.height())
        return;

    fixed_8_120_t ratio = (double)r.height() / r.width();
    fixed_8_120_t ysize = (m_xmax - m_xmin) * (ratio >> 1);
    m_ymin = ((m_ymax + m_ymin) >> 1) - ysize;
    m_ymax = m_ymin + (ysize << 1);
    // DebugPrint(L"SetAspectRatio: ysize=%lf, m_ymin=%lf, m_ymax=%lf\n", (double)ysize, (double)m_ymin, (double)m_ymax);
}

void QMandelbrotWidget::OnZoomChange(const QPoint& point, double zoomMultiplier)
{
    QSize s = size();
    static fixed_8_120_t one = 1;

    // DebugPrint(L"OnZoomChange: Zoom level: %0.10lf\n", (double)m_zoom);
    //fix y coords
    fixed_8_120_t alpha = (double)(point.y()) / (s.height() - 1);
    fixed_8_120_t quarter = (m_ymax - m_ymin) * (1.0 / (zoomMultiplier * 2.0));
    fixed_8_120_t center = alpha * m_ymax + (one - alpha) * m_ymin;
    // DebugPrint(L"OnZoomChange: Y calc: alpha=%0.10lf, quarter=%0.10lf, center=%0.10lf\n", (double)alpha, (double)quarter, (double)center);

    m_ymin = center - quarter;
    m_ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x()) / (s.width() - 1);
    quarter = (m_xmax - m_xmin) * (1.0 / (zoomMultiplier * 2.0));
    center = alpha * m_xmax + (one - alpha) * m_xmin;
    m_xmin = center - quarter;
    m_xmax = center + quarter;

    //DebugPrint(L"OnZoomChange: X calc: alpha=%0.10lf, quarter=%0.10lf, center=%0.10lf\n", (double)alpha, (double)quarter, (double)center);
    //DebugPrint(L"OnZoomChange coords: \n\tX: {%0.10lf, %0.10lf} \n\tY: {%0.10lf, %0.10lf}\n", (double)m_xmin, (double)m_xmax, (double)m_ymin, (double)m_ymax);
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::CreateColorTables()
{
    m_colorTable.clear();
    m_colorTable.resize(m_MaxIter + 1);

    if (m_PaletteType == palGrey) {
        for (int64_t i = 1; i <= m_MaxIter; ++i) {
            int c = 255 - (int)(215.0f * (float)i / (float)m_MaxIter);
            m_colorTable[i] = qRgb(c, c, c);
        }
    }
    else if (m_PaletteType == palVivid) {
        float step = (m_MaxIter <= 256) ? (13.0f / 256.f) : (13.0f / 256.f);
        for (int i = 0; i <= m_MaxIter; ++i) {
            float h = step * i;
            float x = (1.0f - fabs(fmodf(h, 2) - 1.0f));
            switch ((int)floorf(h) % 6) {
            case 0: // 0-60 
                m_colorTable[i] = qRgb(255, int(x * 255), 0);
                break;
            case 1: // 60-120
                m_colorTable[i] = qRgb(int(x * 255), 255, 0);
                break;
            case 2: // 120-180 
                m_colorTable[i] = qRgb(0, 255, int(x * 255));
                break;
            case 3: // 180-240
                m_colorTable[i] = qRgb(0, int(x * 255), 255);
                break;
            case 4: // 240-300
                m_colorTable[i] = qRgb(int(x * 255), 0, 255);
                break;
            case 5: // 300-360
                m_colorTable[i] = qRgb(255, 0, int(x * 255));
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
            m_colorTable[i] = qRgb(b, g, r);
        }
    }

    m_colorTable[0] = qRgb(255, 255, 255);
    m_colorTable[m_MaxIter] = qRgb(0, 0, 0);
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
            m_colorTable[i] = BGR_to_QRgb(val, 0, 255);
            break;
        case 1:
            m_colorTable[i] = BGR_to_QRgb(255, 0, 255 - val);
            break;
        case 2:
            m_colorTable[i] = BGR_to_QRgb(255, val, 0);
            break;
        case 3:
            m_colorTable[i] = BGR_to_QRgb(255 - val, 255, 0);
            break;
        case 4:
            m_colorTable[i] = BGR_to_QRgb(0, 255, val);
            break;
        case 5:
            m_colorTable[i] = BGR_to_QRgb(0, 255 - val, 255);
            break;
        default:
            m_colorTable[i] = qRgb(255, 255, 255);
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
                scanLine[k] = m_colorTable[0];
                continue;
            }
            if (index == (unsigned int)(m_MaxIter - 1)) {
                scanLine[k] = m_colorTable[index];
            }
            else {
                QRgb c1 = m_colorTable[index];
                QRgb c2 = m_colorTable[index + 1];
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

void QMandelbrotWidget::DrawImageDouble(float* pIterations, int64_t w, int64_t h, double x0, double dx, double y0, double dy, double cr, double ci)
{
    const float radius_sq = 2.0F * 2.0F;
    const float LOG2 = logf(2.0F);

    double* xTable = new double[w];
    for (int i = 0; i < w; ++i) xTable[i] = x0 + (double)i * dx;

    bool isJulia = cr != 0.0 || ci != 0.0;

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

            while (modulus < radius_sq && ++iter < m_MaxIter) {
                double tmp = usq - vsq + xc;
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
    const fixed_8_120_t& y0, const fixed_8_120_t& dy, const fixed_8_120_t& cr, const fixed_8_120_t& ci)
{
    const fixed_8_120_t radius_sq = 2 * 2;
    const float LOG2 = logf(2.0F);

    fixed_8_120_t* xTable = new fixed_8_120_t[width];
    for (int i = 0; i < width; ++i) xTable[i] = x0 + dx * i;

    bool isJulia = cr || ci;

#pragma omp parallel for schedule(dynamic)
    for (int l = 0; l < height; ++l) {
        fixed_8_120_t y = y0 + (dy * l);
        fixed_8_120_t usq, vsq, u, v, x, tmp, uv, modulus;
        fixed_8_120_t xc = (isJulia) ? cr : fixed_8_120_t();
        fixed_8_120_t yc = (isJulia) ? ci : y;

        float* pbuff = pIterations + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            x = xTable[k];
            if (isJulia) {
                u = x; v = y; usq = u * u; vsq = v * v; modulus = usq + vsq;
            }
            else {
                u = 0u; v = 0u; usq = 0u; vsq = 0u; xc = x; modulus = 0u;
            }

            while (modulus < radius_sq && ++iter < m_MaxIter) {
                tmp = usq - vsq + xc;
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
    if (m_image.isNull() || m_image.size() != size()) {
        m_image = QImage(size(), QImage::Format_RGB32);
        m_image.fill(Qt::white);

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
        fixed_8_120_t dx = (m_xmax - m_xmin) * (1.0 / w);
        fixed_8_120_t dy = dx;

        fixed_8_120_t cr = 0, ci = 0;
        if (m_SetType == stJulia) {
            cr = m_JuliaCr, ci = m_JuliaCi;
        }

        if (m_precision == Precision::Double || (m_precision == Precision::Auto && m_zoom <= (1ull<<44))) {
            DrawImageDouble(m_Iterations, w, h, (double)m_xmin, (double)dx, (double)m_ymin, (double)dy, cr, ci);
        }
        else {
            DrawImageFixedPoint128(m_Iterations, w, h, m_xmin, dx, m_ymin, dy, cr, ci);
        }

        if (m_PaletteType == palHistogram) {
            CreateHistogram(m_Iterations, w, h);
            CreateColorTableFromHistogram(m_HsvOffset);
        }

        m_NeedToRecompute = false;
    }

    CreateDibFromIterations(m_image, m_Iterations, w, h);
    p.drawImage(rect().topLeft(), m_image);

    // prepare and emit frame stats
    FrameStats stats;
    qint64 elapsed = _paintTimer.elapsed();
    stats.render_time_ms = static_cast<uint32_t>(elapsed > UINT32_MAX ? UINT32_MAX : elapsed);
    stats.zoom = static_cast<float>(m_zoom);
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
        m_zoom *= zoomMultiplier;
        OnZoomChange(event->pos(), zoomMultiplier);
    }
    else if (event->button() == Qt::RightButton) {
        double zoomMultiplier = 0.5;
        if (event->modifiers() & Qt::ControlModifier) {
            zoomMultiplier = (event->modifiers() & Qt::ShiftModifier) ? 0.125 : 0.25;
        }
        m_zoom *= zoomMultiplier;
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
    double logZoom = max(log2(m_zoom), 0);
    
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
    m_JuliaCr = c.real();
    m_JuliaCi = c.imag();
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::setSetType(set_type_t type)
{
    if (type >= stCount || type == m_SetType)
        return;

    m_SetType = type; 
    m_NeedToRecompute = true; 
    update();
}

void QMandelbrotWidget::resetView()
{
    SetDefaultValues();
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::setAnimatePalette(bool animate)
{
    m_animate = animate;
    if (m_animate) {
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
    m_precision = p;
    m_NeedToRecompute = true;
    update();
}

void QMandelbrotWidget::animationTick()
{
    // roll the m_colorTable values

    if (m_PaletteType == palHistogram) {
        m_HsvOffset += 1.0f / 30;
        CreateColorTableFromHistogram(m_HsvOffset);
    }
    else {
        QVector<QRgb> temp = m_colorTable;
        temp[0] = m_colorTable[0];
        for (int i = 1; i < m_MaxIter; ++i) {
            temp[i] = m_colorTable[i + 1];
        }
        temp[m_MaxIter] = m_colorTable[1];

        m_colorTable = std::move(temp);
    }
    repaint();
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
