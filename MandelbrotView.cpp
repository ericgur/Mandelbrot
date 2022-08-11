/***********************************************************************************
7    MIT License

    Copyright (c) 2022 Eric Gur (ericgur@iname.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
************************************************************************************/

// MandelbrotView.cpp : implementation of the CMandelbrotView class
//

#include "pch.h"
#include "Mandelbrot.h"
#include "MandelbrotDoc.h"
#include "MandelbrotView.h"
#include "ComplexSelectDlg.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// debug helper macros
//#define DISABLE_OMP 1
//#define FLOAT128_DEBUG 1
//#define INITIAL_POINT 1
//#define PROFILING 1

#ifdef FLOAT128_DEBUG
#define MAX_ZOOM 0   // for debug purposes
#else
#define MAX_ZOOM (1ull<<44)
#endif

#ifdef INITIAL_POINT
double x_init[] = {-1.1402282717, -1.1400756838};
double y_init[] = {0.2097799200, 0.2098570171};
double zoom_init = 32768;
#endif

//inline COLORREF blendAlpha(COLORREF a, COLORREF b, DWORD alpha)
//{
//    //res = A + (B–A) * alpha
//    COLORREF c1 = (a & 0xFF)     + ((alpha * ((b & 0xFF) - (a & 0xFF))) >> 8);
//    COLORREF c2 = 0;//  (a & 0xFF00) + ((alpha * ((b & 0xFF00) - (a & 0xFF00))) >> 8);
//    COLORREF c3 = 0;// (a & 0xFF0000) + ((alpha * ((b & 0xFF0000) - (a & 0xFF0000))) >> 8);
//    COLORREF res = (c1 & 0xFF) | (c2 & 0xFF00) | (c3 & 0xFF0000);
//    return res;
//}

inline COLORREF blendAlpha(COLORREF colora, COLORREF colorb, DWORD alpha)
{
    COLORREF rb1 = (0x100 - alpha) * (colora & 0xFF00FF);
    COLORREF rb2 = alpha * (colorb & 0xFF00FF);
    COLORREF g1 = (0x100 - alpha) * (colora & 0x00FF00);
    COLORREF g2 = alpha * (colorb & 0x00FF00);
    return (((rb1 + rb2) >> 8) & 0xFF00FF) + (((g1 + g2) >> 8) & 0x00FF00);
}

#ifdef DEBUG
void DebugPrint(const TCHAR* fmt, ...)
{
    TCHAR msgBuffer[2048];
    va_list argList;
    va_start(argList, fmt);
    _vstprintf(msgBuffer, 2048, fmt, argList);
    va_end(argList);
    OutputDebugString(msgBuffer);
}
#else
#define DebugPrint(...)
#endif //#ifdef DEBUG

// CMandelbrotView

IMPLEMENT_DYNCREATE(CMandelbrotView, CView)

BEGIN_MESSAGE_MAP(CMandelbrotView, CView)
    ON_WM_LBUTTONDOWN()
    ON_WM_RBUTTONDOWN()
    ON_WM_MBUTTONDOWN()
    ON_WM_CREATE()
    ON_WM_SIZE()
    ON_WM_SIZING()
    ON_WM_TIMER()
    ON_COMMAND_RANGE(ID_ITERATIONS, ID_ITERATIONS_LAST, OnIterationChange)
    ON_COMMAND_RANGE(ID_SAVEIMAGE_1920X1080, ID_SAVEIMAGE_3840X2160, &CMandelbrotView::OnFileSaveImage)
    ON_COMMAND(ID_VIEW_RESETVIEW, &CMandelbrotView::OnResetView)
    ON_COMMAND(ID_VIEW_ANIMATEPALETTE, &CMandelbrotView::OnAnimatePalette)
    ON_COMMAND_RANGE(ID_VIEW_GREYSCALE, ID_VIEW_HISTOGRAMCOLORING, &CMandelbrotView::OnPaletteChange)
    ON_COMMAND_RANGE(ID_SETTYPE_MANDELBROT, ID_SETTYPE_JULIA, &CMandelbrotView::OnSetTypeSelect)
    ON_COMMAND(ID_SETTYPE_CHOOSEJULIACONSTANT, &CMandelbrotView::OnSetTypeChooseJuliaConstant)
    ON_COMMAND(ID_VIEW_SMOOTHCOLORTRANSITION, &CMandelbrotView::OnSmoothColorTransitions)
    ON_WM_ENTERSIZEMOVE()
    ON_WM_EXITSIZEMOVE()
END_MESSAGE_MAP()


/**
 * @brief CMandelbrotView Constructor
*/
CMandelbrotView::CMandelbrotView()
{
    m_MaxIter = 128;
    m_SmoothLevel = true;
    m_ColorTable32 = nullptr;
    m_Histogram = nullptr;
    m_BmpBits = nullptr;
    m_Iterations = nullptr;
    m_BuffLen = 0;
    m_SetType = stMandelbrot;
    m_JuliaCr = 0.285;
    m_JuliaCi = 0.01;
    m_IsResizing = false;
    m_AnimatePalette = false;
    m_TimerID = 0;

    //fill bitmap header
    memset(&m_BmpInfo, 0, sizeof(m_BmpInfo));
    m_BmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    m_BmpInfo.bmiHeader.biPlanes = 1;
    m_BmpInfo.bmiHeader.biBitCount = 32;
    m_BmpInfo.bmiHeader.biCompression = BI_RGB;
    m_BmpInfo.bmiHeader.biXPelsPerMeter = 100;
    m_BmpInfo.bmiHeader.biYPelsPerMeter = 100;

    //init color table
    m_PaletteType = palGradient;
    CreateColorTables();
    m_NeedToRecompute = true;

    LARGE_INTEGER li;
    QueryPerformanceFrequency(&li);
    m_Frequency = double(li.QuadPart);
}


/**
 * @brief CMandelbrotView Destructor
*/
CMandelbrotView::~CMandelbrotView()
{
    delete[] m_ColorTable32;
    if (m_BmpBits)
        free(m_BmpBits);
    if (m_Iterations)
        delete[] m_Iterations;
    if (m_Histogram)
        delete[] m_Histogram;
}


BOOL CMandelbrotView::PreCreateWindow(CREATESTRUCT& cs)
{
    // TODO: Modify the Window class or styles here by modifying
    //  the CREATESTRUCT cs

    return CView::PreCreateWindow(cs);
}


void CMandelbrotView::OnEnterSizeMove()
{ 
    m_IsResizing = true; 
}


void CMandelbrotView::OnExitSizeMove() 
{ 
    m_IsResizing = false; 
    m_NeedToRecompute = true;
    Invalidate(FALSE);
}


void CMandelbrotView::CreateHistogram(const float* pIterations, int width, int height)
{
    if (m_Histogram)
        delete[] m_Histogram;

    m_Histogram = new int[m_MaxIter + 2];
    ZeroMemory(m_Histogram, sizeof(int) * (m_MaxIter + 2));
    float* hues = new float[m_MaxIter + 1];
    ZeroMemory(hues, sizeof(float) * (m_MaxIter + 1));

#pragma omp parallel
    {
        size_t alloc_size = sizeof(int) * (m_MaxIter + 1);
        int* histogram_private = (int*)_malloca(alloc_size);
        ZeroMemory(histogram_private, alloc_size);
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
    }

    // normalize the iterations
    int total = 0;
    for (int i = 0; i < m_MaxIter; ++i) { // not counting the max iteration score (always black)
        total += m_Histogram[i];
    }

    float hue = 0;
    for (int i = 0; i < m_MaxIter; ++i) {
        hue += (float)m_Histogram[i] / total;
        hues[i] = min(hue, 1.f);
    }
    hues[m_MaxIter] = min(hue, 1.f);

    // create HSV to COLORREF table
    for (int i = 0; i <= m_MaxIter; ++i) {
        float h = hues[i] * 6;
        float x = (1.0f - abs(fmodf(h, 2) - 1.0f));
        switch ((int)floorf(h) % 6) {
        case 0: // 0-60 
            m_ColorTable32[i] = RGB(255, int(x * 255), 0);
            break;
        case 1: // 60-120
            m_ColorTable32[i] = RGB(int(x * 255), 255, 0);
            break;
        case 2: // 120-180 
            m_ColorTable32[i] = RGB(0, 255, int(x * 255));
            break;
        case 3: // 180-240
            m_ColorTable32[i] = RGB(0, int(x * 255), 255);
            break;
        case 4: // 240-300
            m_ColorTable32[i] = RGB(int(x * 255), 0, 255);
            break;
        case 5: // 300-360
            m_ColorTable32[i] = RGB(255, 0, int(x * 255));
            break;
        default: //bug
            m_ColorTable32[i] = RGB(255, 255, 255);
        }
    }

    delete[] hues;
}


void CMandelbrotView::CreateDibFromIterations(COLORREF* pBits, const float* pIterations, int width, int height)
{

    LARGE_INTEGER time_start, time_end;
    QueryPerformanceCounter(&time_start);

#ifndef DISABLE_OMP
#pragma omp parallel for
#endif
    for (int l = 0; l < height; ++l) {
        //point to start of buffer
        COLORREF* pDibPixel = pBits + width * l;
        const float* pIter = pIterations + width * l;

        for (int k = 0; k < width; ++k) {
            float mu = *pIter++;
            if (mu >= m_MaxIter) {
                *(pDibPixel++) = 0;
                continue;
            }

            DWORD index = (DWORD)floor(mu);
            if (index < 0) {
                *(pDibPixel++) = m_ColorTable32[0];
            }
            else {
                COLORREF c1 = m_ColorTable32[index];
                COLORREF c2 = m_ColorTable32[index + 1];
                DWORD alpha = (DWORD)(255.0 * (mu - index));
                if (alpha > 255)
                    alpha = 255;
                *(pDibPixel++) = blendAlpha(c1, c2, alpha);
            }
        }
    }


    QueryPerformanceCounter(&time_end);

#ifdef PROFILING
    DWORD totalTime = DWORD(1000.0 * (time_end.QuadPart - time_start.QuadPart) / m_Frequency);

    CString text;
    text.Format(L"CreateDibFromIterations: (%ims)\n", totalTime);
    OutputDebugString(text);
#endif
}

/**
 * @brief Draw the Mandelbrot image on a DIB surface - uses high precision fixed point
 * @param pIterations: output iterations per pixel
 * @param width: width in pixels
 * @param height: height in pixels
 * @param x0: left most coord
 * @param dx: delta coord between pixels
 * @param y0: top or bottom most coord. Depending if the image is top down or bottom up
 * @param dy: delta coord between pixels, negative for top down DIBs
 * @param cr: Julia constant (Real part)
 * @param ci: Julia constant (Imaginary part)
*/
void CMandelbrotView::DrawImageFixedPoint128(float* pIterations, int width, int height, const fixed_8_120_t& x0, const fixed_8_120_t& dx,
                                             const fixed_8_120_t& y0, const fixed_8_120_t& dy, const fixed_8_120_t& cr, const fixed_8_120_t& ci)
{
    fixed_8_120_t radius = 2, radius_sq = radius * radius;
    const float LOG2 = logf(2.0);

    //create x table
    fixed_8_120_t* xTable = new fixed_8_120_t[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

    bool isJulia = cr || ci ;

#ifndef DISABLE_OMP
#pragma omp parallel for
#endif
    for (int l = 0; l < height; ++l) {
        fixed_8_120_t y = y0 + (dy * l);
        fixed_8_120_t usq, vsq, u, v, x, tmp, uv, modulus;
        fixed_8_120_t xc = (isJulia) ? cr : 0; // no need to do this per pixel for Julia
        fixed_8_120_t yc = (isJulia) ? ci : y;

        //point to start of buffer
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
                u = 0;
                v = 0;
                usq = 0;
                vsq = 0;
                xc = x;
                modulus = 0;
            }

            // complex iterative equation is:
            // C(i) = C(i-1) ^ 2 + C(0)
            while (modulus < radius_sq && ++iter < m_MaxIter) {
                // real
                tmp = usq - vsq + xc;

                // imaginary
                //v = 2.0 * (u * v) + y;
                v = ((u * v) << 1) + yc;
                u = tmp;
                usq = u * u;
                vsq = v * v;
                // check uv vector amplitude is smaller than 2
                modulus = usq + vsq;
            }

            if (m_SmoothLevel && iter < m_MaxIter && iter > 1) {
                *(pbuff++) = (float)(iter + 1) - (logf(logf(sqrtf((float)modulus)))) / LOG2;
            }
            else {
                *(pbuff++) = (float)max(iter, 1);
            }
        }
    }

    delete[] xTable;
}


/**
 * @brief Draw the Mandelbrot or Julia image on a DIB surface - uses double precision floats
 * @param pIterations: output iterations per pixel
 * @param width: width in pixels
 * @param height: height in pixels
 * @param x0: left most coord
 * @param dx: delta coord between pixels
 * @param y0: top or bottom most coord. Depending if the image is top down or bottom up
 * @param dy: delta coord between pixels, negative for top down DIBs
 * @param cr: Julia constant (Real part)
 * @param ci: Julia constant (Imaginary part)
*/
void CMandelbrotView::DrawImageDouble(float* pIterations, int width, int height, double x0, double dx, double y0, double dy, double cr, double ci)
{
    const float radius = 2.0, radius_sq = radius * radius;
    const float LOG2 = logf(2.0);

    //create x table
    double* xTable = new double[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + (double)(i) * dx;
    }

    bool isJulia = cr != 0 || ci != 0;

#ifndef DISABLE_OMP
#pragma omp parallel for
#endif
    for (int l = 0; l < height; ++l) {
        double y = y0 + (dy * l);
        double usq = 0, vsq = 0, u = 0, v = 0;
        double xc = (isJulia) ? cr : 0;
        double yc = (isJulia) ? ci : y;
        double modulus = 0;

        //point to start of buffer
        float* pbuff = pIterations + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            double x = xTable[k];
            // Julia
            if (isJulia) {
                u = x, v = y;
                usq = u * u;
                vsq = v * v;
                modulus = usq + vsq;
            }
            // Mandelbrot
            else {
                u = 0, v = 0; xc = x;
                usq = 0, vsq = 0;
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
                *(pbuff++) = max(mu, 1);
            }
            else {
                *(pbuff++) = (float)max(iter, 1);
            }
        }
    }

    delete[] xTable;
}

/**
 * @brief Callback to redraw the image
 * @param pDC - device context
*/
void CMandelbrotView::OnDraw(CDC* pDC)
{
    CString title;
    CRect rect;
    if (m_ymin == m_ymax) {
        SetAspectRatio();
    }

    GetClientRect(rect);
    const int width = rect.Width(), height = rect.Height();

    if (m_IsResizing) {
        if (m_BmpBits != nullptr) {
            int w = m_BmpInfo.bmiHeader.biWidth;
            int h = m_BmpInfo.bmiHeader.biHeight;
            CBitmap bmp;
            bmp.CreateDiscardableBitmap(pDC, w, h);

            CDC memDC;
            memDC.CreateCompatibleDC(pDC);
            memDC.SelectObject(bmp);
            SetDIBitsToDevice((HDC)(memDC), 0, 0, w, h, 0, 0, 0, h, m_BmpBits, &m_BmpInfo, DIB_RGB_COLORS);

            pDC->StretchBlt(0, 0, width, height, &memDC, 0, 0, w, h, SRCCOPY);
            return;
        }
    }

    // deallocate on size change
    if (m_BmpInfo.bmiHeader.biHeight != height || m_BmpInfo.bmiHeader.biWidth != width) {
        free(m_BmpBits);
        m_BmpBits = nullptr;
        if (m_Iterations != nullptr)
            delete[] m_Iterations;
        m_Iterations = nullptr;
    }

    // allocate new bitmap if needed
    if (nullptr == m_BmpBits) {
        m_BmpInfo.bmiHeader.biHeight = height;
        m_BmpInfo.bmiHeader.biWidth = width;
        m_BuffLen = (size_t)height * width;
        m_BmpBits = (COLORREF*)malloc(m_BuffLen * sizeof(COLORREF));
        m_Iterations = new float[width * height];
        m_NeedToRecompute = true;
    }

    fixed_8_120_t dx = (m_xmax - m_xmin) * (1.0 / width), dy = dx;

    if (m_NeedToRecompute) {
        LARGE_INTEGER time_start, time_end;
        QueryPerformanceCounter(&time_start);

        fixed_8_120_t cr = 0, ci = 0;
        if (m_SetType == stJulia)             {
            cr = m_JuliaCr, ci = m_JuliaCi;
        }
        if (m_zoom > MAX_ZOOM) {
            DrawImageFixedPoint128(m_Iterations, width, height, m_xmin, dx, m_ymin, dy, cr, ci);
        }
        else {
            DrawImageDouble(m_Iterations, width, height, m_xmin, dx, m_ymin, dy, cr, ci);
        }

        if (palHistogram == m_PaletteType) {
            CreateHistogram(m_Iterations, width, height);
        }

        //all done
        QueryPerformanceCounter(&time_end);

        DWORD totalTime = DWORD(1000.0 * (time_end.QuadPart - time_start.QuadPart) / m_Frequency);
        if (m_zoom < 1.0)
            title.Format(L"Zoom x%0.5f (%ims)", m_zoom, totalTime);
        else if (m_zoom > (1 << 16))
            title.Format(L"Zoom x2^%0.0lf (%ims)", log2(m_zoom), totalTime);
        else
            title.Format(L"Zoom x%0.0lf (%ims)", m_zoom, totalTime);

        title += (m_zoom > MAX_ZOOM) ? L" using high precision (slow)" : L" using double precision";
        ((CFrameWnd*)AfxGetMainWnd())->SetWindowText(title);
        m_NeedToRecompute = false;
    }

    // convert iterations to DIB image
    CreateDibFromIterations(m_BmpBits, m_Iterations, width, height);

    ASSERT(m_BmpBits != NULL);
    SetDIBitsToDevice((HDC)(*pDC), 0, 0, width, height, 0, 0, 0, height, m_BmpBits, &m_BmpInfo, DIB_RGB_COLORS);
}

/**
 * @brief Sets the default coordinates of the initial view
*/
void CMandelbrotView::SetDefaultValues()
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


/**
 * @brief Modifies the Y coordinates based on the X coordinates and the resolution
*/
void CMandelbrotView::SetAspectRatio()
{
    CRect rect;
    GetClientRect(rect);
    //use m_xmin, m_xmax and m_rect to determine m_ymin and m_ymax
    //check if window created
    if (0 == rect.Height())
        return;

    fixed_8_120_t ratio = (double)rect.Height() / rect.Width();
    fixed_8_120_t ysize = (m_xmax - m_xmin) * (ratio >> 1);
    m_ymin = ((m_ymax + m_ymin) >> 1) - ysize;
    m_ymax = m_ymin + (ysize << 1);
    DebugPrint(L"SetAspectRatio: ysize=%lf, m_ymin=%lf, m_ymax=%lf\n", (double)ysize, (double)m_ymin, (double)m_ymax);
}


// CMandelbrotView diagnostics

#ifdef _DEBUG
void CMandelbrotView::AssertValid() const
{
    CView::AssertValid();
}


void CMandelbrotView::Dump(CDumpContext& dc) const
{
    CView::Dump(dc);
}


CMandelbrotDoc* CMandelbrotView::GetDocument() const // non-debug version is inline
{
    ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CMandelbrotDoc)));
    return (CMandelbrotDoc*)m_pDocument;
}
#endif //_DEBUG

void CMandelbrotView::OnZoomChange(CPoint& point, double zoomMultiplier)
{
    CRect rect;
    GetClientRect(&rect);
    static fixed_8_120_t one = 1;

    DebugPrint(L"OnZoomChange: Zoom level: %0.10lf\n", (double)m_zoom);
    //fix y coords
    fixed_8_120_t alpha = 1.0 - (double)(point.y) / rect.bottom;
    fixed_8_120_t quarter = (m_ymax - m_ymin) * (1.0 / (zoomMultiplier * 2.0));
    fixed_8_120_t center = alpha * m_ymax + (one - alpha) * m_ymin;
    DebugPrint(L"OnZoomChange: Y calc: alpha=%0.10lf, quarter=%0.10lf, center=%0.10lf\n", (double)alpha, (double)quarter, (double)center);

    m_ymin = center - quarter;
    m_ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (m_xmax - m_xmin) * (1.0 / (zoomMultiplier * 2.0));
    center = alpha * m_xmax + (one - alpha) * m_xmin;
    m_xmin = center - quarter;
    m_xmax = center + quarter;

    DebugPrint(L"OnZoomChange: X calc: alpha=%0.10lf, quarter=%0.10lf, center=%0.10lf\n", (double)alpha, (double)quarter, (double)center);
    DebugPrint(L"OnZoomChange coords: \n\tX: {%0.10lf, %0.10lf} \n\tY: {%0.10lf, %0.10lf}\n", (double)m_xmin, (double)m_xmax, (double)m_ymin, (double)m_ymax);
    SetAspectRatio();
    m_NeedToRecompute = true;
    Invalidate(FALSE);
}

// CMandelbrotView message handlers


/**
 * @brief Called on Mouse Left click. Zooms in the image by 2x or 4x (if CTRL key is pressed)
 * @param nFlags holds which mosue or special keys where pressed befor this event
 * @param point Mouse coord within the client window
*/
void CMandelbrotView::OnLButtonDown(UINT nFlags, CPoint point)
{
    double zoomMultiplier = 2.0;
    if (nFlags & MK_CONTROL) {
        zoomMultiplier = (nFlags & MK_SHIFT) ? 8.0 : 4.0;
    }

    m_zoom *= zoomMultiplier;
    OnZoomChange(point, zoomMultiplier);
}


/**
 * @brief Called on Mouse Right click. Zooms out the image by 2x or 4x (if CTRL key is pressed)
 * @param nFlags holds which mosue or special keys where pressed befor this event
 * @param point Mouse coord within the client window
*/
void CMandelbrotView::OnRButtonDown(UINT nFlags, CPoint point)
{
    double zoomMultiplier = 0.5;
    if (nFlags & MK_CONTROL) {
        zoomMultiplier = (nFlags & MK_SHIFT) ? 0.125 : 0.25;
    }

    m_zoom *= zoomMultiplier;
    OnZoomChange(point, zoomMultiplier);
}


/**
 * @brief Called on Mouse Middle click. Resets the view.
 * @param nFlags holds which mosue or special keys where pressed befor this event (not used)
 * @param point Mouse coord within the client window (not used)
*/

void CMandelbrotView::OnMButtonDown(UINT nFlags, CPoint point)
{
    UNREFERENCED_PARAMETER(nFlags);
    UNREFERENCED_PARAMETER(point);
    OnResetView();
}


/**
 * @brief Reset the view to default coordinates and zoom
*/
void CMandelbrotView::OnResetView()
{
    SetDefaultValues();
    m_NeedToRecompute = true;
    Invalidate(FALSE);
}

/**
 * @brief Creates iteration to color mapping depending on the max iteration count and user selected palette type
 *        Histogram coloring is calculated on the fly
*/
void CMandelbrotView::CreateColorTables()
{
    if (m_ColorTable32)
        delete[] m_ColorTable32;

    m_ColorTable32 = new COLORREF[m_MaxIter + 2];
    if (m_PaletteType == palGrey) {
        for (size_t i = 1; i <= m_MaxIter; ++i) {
            int c = 255 - (int)(215.0f * (float)i / (float)m_MaxIter);
            m_ColorTable32[i] = RGB(c, c, c);
        }
    }
    else if (m_PaletteType == palVivid) {
        float step = 6 * (m_MaxIter <= 256) ? (13.0f / 256.f) : (13.0f / 256.f);
        for (int i = 0; i <= m_MaxIter; ++i) {
            float h = step * i;
            float x = (1.0f - abs(fmodf(h, 2) - 1.0f));
            switch ((int)floorf(h) % 6) {
            case 0: // 0-60 
                m_ColorTable32[i] = RGB(255, int(x * 255), 0);
                break;
            case 1: // 60-120
                m_ColorTable32[i] = RGB(int(x * 255), 255, 0);
                break;
            case 2: // 120-180 
                m_ColorTable32[i] = RGB(0, 255, int(x * 255));
                break;
            case 3: // 180-240
                m_ColorTable32[i] = RGB(0, int(x * 255), 255);
                break;
            case 4: // 240-300
                m_ColorTable32[i] = RGB(int(x * 255), 0, 255);
                break;
            case 5: // 300-360
                m_ColorTable32[i] = RGB(255, 0, int(x * 255));
                break;
            }
        }
    }
    else if (m_PaletteType == palGradient) {
        UINT r = 0;
        UINT g = 20;
        UINT b = 255;

        for (size_t i = 1; i <= m_MaxIter; ++i) {
            r = (r + 3) & 0xFF;
            g = (g + 5) & 0xFF;
            b = (b - 3) & 0xFF;
            m_ColorTable32[i] = RGB(b, g, r);
        }
    }

    m_ColorTable32[0] = RGB(255, 255, 255);
}


int CMandelbrotView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
    if (CView::OnCreate(lpCreateStruct) == -1)
        return -1;

    SetDefaultValues();
    Invalidate(FALSE);
    return 0;
}

/**
 * @brief - callback when user selects Mandelbrot or Julia sets
 * @param nID - resource ID of the menu item
*/
void CMandelbrotView::OnSetTypeSelect(UINT nID)
{
    set_type_t set_type = stMandelbrot;
    switch (nID) {
    case ID_SETTYPE_MANDELBROT:
        set_type = stMandelbrot;
        break;
    case ID_SETTYPE_JULIA:
        set_type = stJulia;
        break;
    default:
        DebugPrint(L"OnSetTypeSelect: Invalid nID received %d\n", nID);
        break;
    }

    // no change
    if (set_type == m_SetType)
        return;
    
    m_SetType = set_type;
    CMenu* menu = AfxGetMainWnd()->GetMenu();
    menu->CheckMenuRadioItem(ID_SETTYPE_MANDELBROT, ID_SETTYPE_JULIA, nID, MF_BYCOMMAND);
    
    // Reset coordinates
    SetDefaultValues();
    m_NeedToRecompute = true;
    Invalidate(FALSE);
}


/**
 * @brief Callback for menu->iteration selection change. Changes the iteration count and initiates a redraw of the image
 * @param nID Resource ID of the menu item selected
*/
void CMandelbrotView::OnIterationChange(UINT nID)
{
    if (nID < ID_ITERATIONS || nID > ID_ITERATIONS_LAST) {
        CString str;
        str.Format(L"OnIterationChange: Invalid nID received as argument: %d\n", nID);
        OutputDebugString(str);
    }

    CMenu* menu = AfxGetMainWnd()->GetMenu();
    CString value;
    menu->CheckMenuRadioItem(ID_ITERATIONS, ID_ITERATIONS_LAST, nID, MF_BYCOMMAND);
    menu->GetMenuString(nID, value, MF_BYCOMMAND);
    m_MaxIter = _ttoi(value);
    CreateColorTables();
    m_NeedToRecompute = true;
    Invalidate(FALSE);
}


void CMandelbrotView::OnPaletteChange(UINT nID)
{
    if (nID < ID_VIEW_GREYSCALE || nID > ID_VIEW_HISTOGRAMCOLORING) {
        CString str;
        str.Format(L"OnPaletteChange: Invalid nID received as argument: %d\n", nID);
        OutputDebugString(str);
    }

    CMenu* menu = AfxGetMainWnd()->GetMenu();
    CString value;
    menu->CheckMenuRadioItem(ID_VIEW_GREYSCALE, ID_VIEW_HISTOGRAMCOLORING, nID, MF_BYCOMMAND);
    switch (nID) {
    case ID_VIEW_GREYSCALE:
        m_PaletteType = palGrey; break;
    case ID_VIEW_VIVIDCOLORS:
        m_PaletteType = palVivid; break;
    case ID_VIEW_GRADIENTS:
        m_PaletteType = palGradient; break;
    case ID_VIEW_HISTOGRAMCOLORING:
        m_PaletteType = palHistogram; break;
    }

    CreateColorTables();
    m_NeedToRecompute = (m_PaletteType == palHistogram);
    Invalidate(FALSE);
}


/**
 * @brief Callback for ID_FILE_SAVE_IMAGE menu item. Save the image to disk. Uses Same X coordinates but modifies the Y coordinates to fit the resolution.
*/
void CMandelbrotView::OnFileSaveImage(UINT nID)
{
    int width = 0, height = 0; // TODO: add more resolutions.
    switch (nID) {
    case ID_SAVEIMAGE_1920X1080:
        width = 1920, height = 1080;
        break;
    case ID_SAVEIMAGE_2560X1440:
        width = 2560, height = 1440;
        break;
    case ID_SAVEIMAGE_3840X2160:
        width = 3840, height = 2160;
        break;
    default:
        DebugPrint(L"OnSetTypeSelect: Invalid nID received %d\n", nID);
        break;
    }

    CImage image;
    image.Create(width, -height, 32);
    CFileDialog dlg(FALSE,                         //bOpenFileDialog,
                    L"png",                        //LPCTSTR lpszDefExt = NULL,
                    L"untitled",                   //LPCTSTR lpszFileName = NULL,
                    OFN_HIDEREADONLY | OFN_ENABLESIZING
                    | OFN_OVERWRITEPROMPT,         //DWORD dwFlags = OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT,
                    L"PNG files (*.png)|*.png||",  //LPCTSTR lpszFilter = NULL,
                    this);                         //CWnd* pParentWnd = NULL,

    if (IDCANCEL == dlg.DoModal())
        return; //user pressed cancel

    CString filename = dlg.GetPathName();
    if (filename.IsEmpty())
        return;


    COLORREF* pBits = (COLORREF*)image.GetPixelAddress(0, 0);
    float* pIterations = new float[width * height];
    fixed_8_120_t dx = (m_xmax - m_xmin) * (1.0 / width);
    fixed_8_120_t ratio = (double)(height) / width;
    fixed_8_120_t ysize = (m_xmax - m_xmin) * (ratio >> 1);
    fixed_8_120_t ymax = ((m_ymax + m_ymin) >> 1) + ysize;
    fixed_8_120_t cr = 0, ci = 0;
    if (m_SetType == stJulia) {
        cr = m_JuliaCr, ci = m_JuliaCi;
    }

    if (m_zoom > MAX_ZOOM) {
        DrawImageFixedPoint128(pIterations, width, height, m_xmin, dx, ymax, -dx, cr, ci);
    }
    else {
        DrawImageDouble(pIterations, width, height, m_xmin, dx, ymax, -dx, cr, ci);
    }

    CreateDibFromIterations(pBits, pIterations, width, height);

    HRESULT hr = image.Save(filename, Gdiplus::ImageFormatPNG);
    if (!SUCCEEDED(hr)) { 
        AfxMessageBox(L"Failed to save image", MB_ICONINFORMATION);
    }

    delete[] pIterations;
}

/**
 * @brief Callback for ID_SETTYPE_CHOOSEJULIACONSTANT, triggers a dialog that allows the user to select a new Julia set constant.
*/
void CMandelbrotView::OnSetTypeChooseJuliaConstant()
{
    ComplexSelectDlg dlg(this);
    dlg.real.Format(L"%lf", m_JuliaCr);
    dlg.imag.Format(L"%lf", m_JuliaCi);
    
    if (IDCANCEL == dlg.DoModal())
        return; //user pressed cancel
    
    m_JuliaCr = _ttof(dlg.real);
    m_JuliaCi = _ttof(dlg.imag);
    m_NeedToRecompute = true;
    Invalidate(FALSE);
}


/**
 * @brief Callback for ID_VIEW_ANIMATEPALETTE, starts/stops palette animation
*/
void CMandelbrotView::OnAnimatePalette()
{
    CMenu* menu = AfxGetMainWnd()->GetMenu();
    UINT state = menu->GetMenuState(ID_VIEW_ANIMATEPALETTE, MF_BYCOMMAND);

    if (state & MF_CHECKED) {
        menu->CheckMenuItem(ID_VIEW_ANIMATEPALETTE, MF_UNCHECKED | MF_BYCOMMAND);
        m_AnimatePalette = false;
        if (m_PaletteType == palHistogram) {
            m_NeedToRecompute = true;
        }
        else {
            CreateColorTables(); // reset color tables
        }

        if (m_TimerID != 0) {
            KillTimer(m_TimerID);
            m_TimerID = 0;
        }
    }
    else {
        menu->CheckMenuItem(ID_VIEW_ANIMATEPALETTE, MF_CHECKED | MF_BYCOMMAND);
        m_AnimatePalette = true;
        m_TimerID = SetTimer(1, 70, NULL); // timer causes WM_TIMER message
    }

    Invalidate(FALSE);
}

/**
 * @brief Callback for WM_TIMER events.
 * @param nIDEvent Timer ID
*/
void CMandelbrotView::OnTimer(UINT_PTR nIDEvent)
{
    // Call base class
    //CView::OnTimer(nIDEvent);

    // roll the m_ColorTable32 values
    if (nIDEvent == m_TimerID) {
        COLORREF* pBuff = new COLORREF[m_MaxIter + 1];
        for (int i = 1; i < m_MaxIter; ++i) {
            pBuff[i] = m_ColorTable32[i + 1];
        }
        pBuff[0] = m_ColorTable32[0];
        pBuff[m_MaxIter] = m_ColorTable32[1];

        // Note - since using WM_TIMER message (always on the main thread), there's no need to protect with a mutex.
        delete[] m_ColorTable32;
        m_ColorTable32 = pBuff;
        Invalidate(FALSE);
    }
}

/**
 * @brief Callback for ID_VIEW_SMOOTHCOLORTRANSITION, select linear interpolation of colors.
*/
void CMandelbrotView::OnSmoothColorTransitions()
{
    CMenu* menu = AfxGetMainWnd()->GetMenu();
    UINT state = menu->GetMenuState(ID_VIEW_SMOOTHCOLORTRANSITION, MF_BYCOMMAND);

    if (state & MF_CHECKED) {
        menu->CheckMenuItem(ID_VIEW_SMOOTHCOLORTRANSITION, MF_UNCHECKED | MF_BYCOMMAND);
        m_SmoothLevel = false;
    }
    else {
        menu->CheckMenuItem(ID_VIEW_SMOOTHCOLORTRANSITION, MF_CHECKED | MF_BYCOMMAND);
        m_SmoothLevel = true;
    }

    m_NeedToRecompute = true;
    Invalidate(FALSE);
}
