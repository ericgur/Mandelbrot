#pragma once

#include <complex>
#include <QWidget>
#include <QChronoTimer>

// include fixed point implementation from project root
#include "fixed_point128.h"
using namespace fp128;

typedef fixed_point128<8> fixed_8_120_t;
struct FrameStats {
    uint32_t render_time_ms;
    float zoom;
    QSize size;
    int32_t max_iterations;
};

class QMandelbrotWidget : public QWidget
{
    Q_OBJECT
    Q_DISABLE_COPY_MOVE(QMandelbrotWidget)
public:
    static inline constexpr int64_t max_iterations = 4096;
    static inline constexpr int64_t min_iterations = 128;
    static inline constexpr double logMaxZoom = 113.0; // max of fixed point 128 bit
    enum class Precision { Auto, Double, FixedPoint128 };
    enum palette_t { palGrey, palGradient, palVivid, palHistogram };
    enum set_type_t { stMandelbrot, stJulia, stCount };

    explicit QMandelbrotWidget(QWidget* parent = nullptr);
    virtual ~QMandelbrotWidget();
    void saveImage(int width, int height);
    std::complex<double> juliaConstant() const { return m_JuliaConstant; }
    void setSetType(set_type_t type);

signals:
    void renderDone(FrameStats stats);

public slots:
    void setJuliaConstant(const std::complex<double>& c);
    void resetView();
    void zoomIn();
    void zoomOut();
    void setAnimatePalette(bool animate);
    void setPrecision(Precision p);
    void setMaximumIterations(int64_t maxIter);
    void animationTick();
    void panHorizontal(double amount);
    void panVertical(double amount);

protected:
    void paintEvent(QPaintEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;

private:
    int64_t calcAutoIterationLimits();

    // View state (ported from CMandelbrotView)
    fixed_8_120_t m_Xmin, m_Xmax, m_Ymin, m_Ymax;
    double m_ZoomLevel;
    double m_ZoomIncrement = 2.0;
    int64_t m_MaxIter = 128;
    bool m_AutoIterations = false;

    // image and buffers
    QImage m_ImageCache;
    float* m_Iterations = nullptr; // will allocate when needed
    bool m_NeedToRecompute = true;

    // color / palette data
    QVector<QRgb> m_ColorTable;
    int* m_Histogram = nullptr;
    float m_HsvOffset = 0;
    palette_t m_PaletteType = palGradient;
    bool m_SmoothLevel = true;

    // UI flags
    Precision m_Precision = Precision::Auto;
    bool m_Animate = false;

    // set type and julia constants
    set_type_t m_SetType = stMandelbrot;
    std::complex<double> m_JuliaConstant{ 0.285, 0.01 };
    
    // Timer for animation
    QChronoTimer m_Timer;

    // helpers
    void SetDefaultValues();
    void SetAspectRatio();
    void OnZoomChange(const QPoint& point, double zoomMultiplier);

    // rendering helpers (ported)
    void CreateColorTables();
    void CreateColorTableFromHistogram(float offset);
    void CreateHistogram(const float* pIterations, int64_t width, int64_t height);
    void CreateDibFromIterations(QImage& img, const float* pIterations, int64_t width, int64_t height);
    void DrawImageDouble(float* pIterations, int64_t width, int64_t height, double x0, double dx, double y0, double dy);
    void DrawImageFixedPoint128(float* pIterations, int64_t width, int64_t height, const fixed_8_120_t& x0, const fixed_8_120_t& dx, const fixed_8_120_t& y0,
                       const fixed_8_120_t& dy);
};
