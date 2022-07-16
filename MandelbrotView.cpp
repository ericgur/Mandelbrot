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

//#define DEEP_DEBUG 1

#ifdef DEEP_DEBUG
#define MAX_ZOOM (1ull<<1)   // for debug purposes
#else
#define MAX_ZOOM (1ull<<44)
#endif


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
    ON_COMMAND(ID_VIEW_GREYSCALE, OnGreyScale)
    ON_COMMAND_RANGE(ID_ITERATIONS, ID_ITERATIONS_LAST, OnIterationChange)
    ON_COMMAND(ID_FILE_SAVE_IMAGE, &CMandelbrotView::OnFileSaveImage)
    ON_COMMAND_RANGE(ID_SETTYPE_MANDELBROT, ID_SETTYPE_JULIA, &CMandelbrotView::OnSetTypeSelect)
    ON_COMMAND(ID_SETTYPE_CHOOSEJULIACONSTANT, &CMandelbrotView::OnSetTypeChooseJuliaConstant)
END_MESSAGE_MAP()


/**
 * @brief CMandelbrotView Constructor
*/
CMandelbrotView::CMandelbrotView()
{
    m_MaxIter = 128;
    m_SmoothLevel = true;

    m_ColorTable32 = NULL;
    m_BmpBits = NULL;
    m_BuffLen = 0;
    m_SetType = stMandelbrot;
    m_JuliaCr = 0.285;
    m_JuliaCi = 0.01;
    // Other interesting values :
    // c = complex(-0.7269, 0.1889)
    // c = complex(-0.8, 0.156)
    // c = complex(-0.4, 0.6)

    // set default precision for MPIR
    // TODO: change this depending on zoom level. 64 is good for ~x60
    #ifdef USE_MPIR
    mpf_set_default_prec(64);
    #endif

    //fill bitmap header
    memset(&m_BmpInfo, 0, sizeof(m_BmpInfo));
    m_BmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    m_BmpInfo.bmiHeader.biPlanes = 1;
    m_BmpInfo.bmiHeader.biBitCount = 32;
    m_BmpInfo.bmiHeader.biCompression = BI_RGB;
    m_BmpInfo.bmiHeader.biXPelsPerMeter = 100;
    m_BmpInfo.bmiHeader.biYPelsPerMeter = 100;

    //init color table
    m_GreyScale = false;
    CreateColorTables();
    m_NeedToRedraw = true;

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
}


BOOL CMandelbrotView::PreCreateWindow(CREATESTRUCT& cs)
{
    // TODO: Modify the Window class or styles here by modifying
    //  the CREATESTRUCT cs

    return CView::PreCreateWindow(cs);
}

void CMandelbrotView::DrawImageFixedPoint128(COLORREF* pBits, int width, int height, const fixed_8_120_t& x0, const fixed_8_120_t& dx,
                                             const fixed_8_120_t& y0, const fixed_8_120_t& dy, const fixed_8_120_t& cr, const fixed_8_120_t& ci)
{
    fixed_8_120_t radius = 2.0, radius_sq = radius * radius;
    const double LOG2 = log(2.0);

    //create x table
    fixed_8_120_t* xTable = new fixed_8_120_t[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

    bool isJulia = cr || ci ;

#ifndef DEEP_DEBUG
#pragma omp parallel for
#endif
    for (int l = 0; l < height; ++l) {
        fixed_8_120_t y = y0 + (dy * l);
        fixed_8_120_t usq = 0.0, vsq = 0.0, u = 0.0, v = 0.0, x, tmp = 0.0, uv, modulus;
        fixed_8_120_t xc = (isJulia) ? cr : 0.0; // no need to do this per pixel
        fixed_8_120_t yc = (isJulia) ? ci : y;

        //point to start of buffer
        COLORREF* pbuff = pBits + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            x = xTable[k];
            if (isJulia) {
                u = x;
                v = y;
            }
            else {
                u = 0.0;
                v = 0.0;
                xc = x;
            }
            usq = u * u;
            vsq = v * v;
            modulus = usq + vsq;

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

            if (m_SmoothLevel && iter < m_MaxIter && iter > 0) {
                double mu = (double)iter + 1 - (log(log(sqrt(modulus)))) / LOG2;
                DWORD index = (DWORD)floor(mu);
                COLORREF c1 = m_ColorTable32[index];
                COLORREF c2 = m_ColorTable32[index + 1];
                DWORD alpha = (DWORD)(255.0 * (mu - index));
                if (alpha > 255)
                    alpha = 255;
                COLORREF color = blendAlpha(c1, c2, alpha);
                *(pbuff++) = color;
            }
            else {
                *(pbuff++) = m_ColorTable32[iter];
            }
        }
    }

    delete[] xTable;
}


/**
 * @brief Draw the Mandelbrot image on a DIB surface - uses high precision floats
 * @param pBits: output DIB surface
 * @param width: width in pixels
 * @param height: height in pixels
 * @param x0: left most coord
 * @param dx: delta coord between pixels
 * @param y0: top or bottom most coord. Depending if the image is top down or bottom up
 * @param dy: delta coord between pixels, negative for top down DIBs
 * @param cr: Julia constant (Real part)
 * @param ci: Julia constant (Imaginary part)
*/
#ifdef USE_MPIR
void CMandelbrotView::DrawImageMPIR(COLORREF* pBits, int width, int height, const mpf_class& x0, const mpf_class& dx,
                                    const mpf_class& y0, const mpf_class& dy, const mpf_class& cr, const mpf_class& ci)
{
    mpf_class radius = 2.0, radius_sq = radius * radius;
    const double LOG2 = log(2.0);

    //create x table
    mpf_class* xTable = new mpf_class[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

    bool isJulia = cr != 0 || ci != 0;

#pragma omp parallel for
    for (int l = 0; l < height; ++l) {
        mpf_class y = y0 + (dy * l);
        mpf_class usq = 0, vsq = 0, u = 0, v = 0, x, tmp = 0, uv, modulus;
        mpf_class xc = (isJulia) ? cr : 0; // no need to do this per pixel
        mpf_class yc = (isJulia) ? ci : y;

        //point to start of buffer
        COLORREF* pbuff = pBits + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            x = xTable[k];
            if (isJulia) {
                u = x;
                v = y;
            }
            else {
                u = 0;
                v = 0;
                xc = x;
            }
            usq = u;
            usq *= u;
            vsq = v;
            vsq *= v;
            modulus = usq;
            modulus += vsq;

            // complex iterative equation is:
            // C(i) = C(i-1) ^ 2 + C(0)
            while (modulus < radius_sq && ++iter < m_MaxIter) {
                // real
                tmp = usq;
                tmp -= vsq;
                tmp += xc;

                // imaginary
                //v = 2.0 * (u * v) + y;
                v *= u;
                v *= 2;
                v += yc;
                u = tmp;
                vsq = v;
                vsq *= v;
                usq = u;
                usq *= u;

                // check uv vector amplitude is smaller than 2
                modulus = vsq;
                modulus += usq;
            }

            if (m_SmoothLevel && iter < m_MaxIter && iter > 0) {
                double mu = (double)iter + 1 - (log(log(sqrt(modulus.get_d())))) / LOG2;
                DWORD index = (DWORD)floor(mu);
                COLORREF c1 = m_ColorTable32[index];
                COLORREF c2 = m_ColorTable32[index + 1];
                DWORD alpha = (DWORD)(255.0 * (mu - index));
                if (alpha > 255)
                    alpha = 255;
                COLORREF color = blendAlpha(c1, c2, alpha);
                *(pbuff++) = color;
            }
            else {
                *(pbuff++) = m_ColorTable32[iter];
            }
        }
    }

    delete[] xTable;
}
#endif //USE_MPIR

/**
 * @brief Draw the Mandelbrot or Julia image on a DIB surface - uses double precision floats
 * @param pBits: output DIB surface
 * @param width: width in pixels
 * @param height: height in pixels
 * @param x0: left most coord
 * @param dx: delta coord between pixels
 * @param y0: top or bottom most coord. Depending if the image is top down or bottom up
 * @param dy: delta coord between pixels, negative for top down DIBs
 * @param cr: Julia constant (Real part)
 * @param ci: Julia constant (Imaginary part)
*/
void CMandelbrotView::DrawImage(COLORREF* pBits, int width, int height, double x0, double dx, double y0, double dy, double cr, double ci)
{
    const double radius = 2.0, radius_sq = radius * radius;
    const double LOG2 = log(2.0);

    //create x table
    double* xTable = new double[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + (double)(i) * dx;
    }

    bool isJulia = cr != 0 || ci != 0;

#pragma omp parallel for
    for (int l = 0; l < height; ++l) {
        double y = y0 + (dy * l);
        double usq = 0, vsq = 0, u = 0, v = 0;
        double xc = (isJulia) ? cr : 0;
        double yc = (isJulia) ? ci : y;

        //point to start of buffer
        COLORREF* pbuff = pBits + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            double x = xTable[k];
            if (isJulia) {
                u = x, v = y;
            }
            else {
                u = 0, v = 0; xc = x;
            }
            usq = u * u;
            vsq = v * v;
            double modulus = usq + vsq;

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
                double mu = (double)iter + 1 - (log(log(sqrt(modulus)))) / LOG2;
                DWORD index = (DWORD)floor(mu);
                COLORREF c1 = m_ColorTable32[index];
                COLORREF c2 = m_ColorTable32[index + 1];
                DWORD alpha = (DWORD)(255.0 * (mu - index));
                if (alpha > 255)
                    alpha = 255;
                COLORREF color = blendAlpha(c1, c2, alpha);
                *(pbuff++) = color;
            }
            else {
                *(pbuff++) = m_ColorTable32[iter];
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
    SetAspectRatio();

    CString title;

    CRect rect;
    GetClientRect(rect);
    const int width = rect.Width(), height = rect.Height();

    //deallocate on size change
    if (m_BmpInfo.bmiHeader.biHeight != height || m_BmpInfo.bmiHeader.biWidth != width) {
        free(m_BmpBits);
        m_BmpBits = NULL;
    }

    //allocate new bitmap if needed
    if (NULL == m_BmpBits) {
        m_BmpInfo.bmiHeader.biHeight = height;
        m_BmpInfo.bmiHeader.biWidth = width;
        m_BuffLen = (size_t)height * width;
        m_BmpBits = (COLORREF*)malloc(m_BuffLen * sizeof(COLORREF));
        m_NeedToRedraw = true;
    }

    fixed_8_120_t dx((m_xmax - m_xmin) * fixed_8_120_t(1.0/ (double)width)), dy = dx;

    if (m_NeedToRedraw) {
        LARGE_INTEGER time_start, time_end;
        QueryPerformanceCounter(&time_start);

        fixed_8_120_t cr = 0.0, ci = 0.0;
        if (m_SetType == stJulia)             {
            cr = m_JuliaCr, ci = m_JuliaCi;
        }
        if (m_zoom > MAX_ZOOM) {
            DrawImageFixedPoint128(m_BmpBits, width, height, m_xmin, dx, m_ymin, dy, cr, ci);
            //DrawImageMPIR(m_BmpBits, width, height, m_xmin, dx, m_ymin, dy, cr, ci);
        }
        else {
            DrawImage(m_BmpBits, width, height, m_xmin, dx, m_ymin, dy, cr, ci);
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
        m_NeedToRedraw = false;
    }

    ASSERT(m_BmpBits != NULL);
    SetDIBitsToDevice((HDC)(*pDC), 0, 0, width, height, 0, 0, 0, height, m_BmpBits, &m_BmpInfo, DIB_RGB_COLORS);
}

/**
 * @brief Sets the default coordinates of the initial view
*/
void CMandelbrotView::SetDefaultValues()
{
    m_xmax = 2.5;
    m_xmin = -m_xmax;
    m_ymax = m_ymin = 0.0;
    m_zoom = 1;
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
    fixed_8_120_t ratio = (double)(rect.Height()) / (double)(rect.Width());
    fixed_8_120_t ysize((m_xmax - m_xmin) * (ratio * fixed_8_120_t(0.5)));
    m_ymin = ((m_ymax + m_ymin) * fixed_8_120_t(0.5)) - ysize;
    m_ymax = m_ymin + (fixed_8_120_t(2.0) * ysize);
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


// CMandelbrotView message handlers

/**
 * @brief Called on Mouse Left click. Zooms in the image by 2x or 4x (if CTRL key is pressed)
 * @param nFlags holds which mosue or special keys where pressed befor this event
 * @param point Mouse coord within the client window
*/
void CMandelbrotView::OnLButtonDown(UINT nFlags, CPoint point)
{
    double zoomMultiplier = 2.0;
    if (nFlags & MK_CONTROL)
        zoomMultiplier = 4.0;

    m_zoom *= zoomMultiplier;
    CRect rect;
    GetClientRect(&rect);

    fixed_8_120_t alpha(1.0 - ((double)(point.y) / (double)rect.bottom));
    
    //fix y coords
    fixed_8_120_t quarter((m_ymax - m_ymin) * fixed_8_120_t(1.0 / (zoomMultiplier * 2.0)));
    fixed_8_120_t center(alpha * m_ymax + (fixed_8_120_t(1.0) - alpha) * m_ymin);
    DebugPrint(L"OnLButtonDown: alpha=%lf, quarter=%lf, center=%lf\n", (double)alpha, (double)quarter, (double)center);

    m_ymin = center - quarter;
    m_ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (m_xmax - m_xmin) * fixed_8_120_t(1.0 / (zoomMultiplier * 2.0));
    center = alpha * m_xmax + (fixed_8_120_t(1.0) - alpha) * m_xmin;
    m_xmin = center - quarter;
    m_xmax = center + quarter;

    m_NeedToRedraw = true;
    Invalidate(FALSE);

    // CView::OnLButtonDown(nFlags, point);
}


/**
 * @brief Called on Mouse Right click. Zooms out the image by 2x or 4x (if CTRL key is pressed)
 * @param nFlags holds which mosue or special keys where pressed befor this event
 * @param point Mouse coord within the client window
*/
void CMandelbrotView::OnRButtonDown(UINT nFlags, CPoint point)
{
    double zoomMultiplier = 2.0;
    if (nFlags & MK_CONTROL)
        zoomMultiplier = 4.0;

    CRect rect;
    GetClientRect(&rect);

    m_zoom /= zoomMultiplier;

    fixed_8_120_t alpha(1.0 - ((double)(point.y) / (double)rect.bottom));
    //fix y coords
    fixed_8_120_t quarter((m_ymax - m_ymin) * fixed_8_120_t(zoomMultiplier / 2.0));
    fixed_8_120_t center(alpha * m_ymax + (fixed_8_120_t(1.0) - alpha) * m_ymin);
    m_ymin = center - quarter;
    m_ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (m_xmax - m_xmin) * fixed_8_120_t(zoomMultiplier / 2.0);
    center = alpha * m_xmax + (fixed_8_120_t(1.0) - alpha) * m_xmin;
    m_xmin = center - quarter;
    m_xmax = center + quarter;

    m_NeedToRedraw = true;
    Invalidate(FALSE);

    //CView::OnRButtonDown(nFlags, point);
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
    SetDefaultValues();
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}


/**
 * @brief Creates iteration to color mapping depending on the max iteration count
*/
void CMandelbrotView::CreateColorTables()
{
    if (m_ColorTable32)
        delete[] m_ColorTable32;

    m_ColorTable32 = new COLORREF[m_MaxIter + 2];
    if (m_GreyScale) {
        for (size_t i = 1; i <= m_MaxIter; ++i) {
            int c = 255 - (int)(255.0f * (float)i / (float)m_MaxIter);
            m_ColorTable32[i] = RGB(c, c, c);
        }
    }
    else {
        for (size_t i = 1; i <= m_MaxIter; ++i) {
            size_t c = m_MaxIter - i;
            m_ColorTable32[i] = RGB(c * 4 & 255, c * 6 & 255, c * 3 & 255);
        }
    }

    m_ColorTable32[0] = RGB(255, 255, 255);
}


int CMandelbrotView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
    if (CView::OnCreate(lpCreateStruct) == -1)
        return -1;

    SetDefaultValues();
    SetAspectRatio();

    return 0;
}


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
    m_NeedToRedraw = true;
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
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}


/**
 * @brief Callback for the ID_VIEW_GREYSCALE command. Toggles between grey and color images.
*/
void CMandelbrotView::OnGreyScale()
{
    CMenu* menu = AfxGetMainWnd()->GetMenu();

    int state = menu->GetMenuState(ID_VIEW_GREYSCALE, MF_BYCOMMAND);

    if (state & MF_CHECKED) {
        menu->CheckMenuItem(ID_VIEW_GREYSCALE, MF_UNCHECKED | MF_BYCOMMAND);
        m_GreyScale = false;
    }
    else {
        menu->CheckMenuItem(ID_VIEW_GREYSCALE, MF_CHECKED | MF_BYCOMMAND);
        m_GreyScale = true;
    }

    CreateColorTables();
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}


/**
 * @brief Callback for ID_FILE_SAVE_IMAGE menu item. Save the image to disk. Uses Same X coordinates but modifies the Y coordinates to fit the resolution.
*/
void CMandelbrotView::OnFileSaveImage()
{
    int width = 2560, height = 1440; // TODO: add more resolutions.

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

    fixed_8_120_t dx = (m_xmax - m_xmin) * fixed_8_120_t(1.0 / width);
    fixed_8_120_t ratio = (double)(height) / (double)(width);
    fixed_8_120_t ysize((m_xmax - m_xmin) * (ratio * fixed_8_120_t(0.5)));
    fixed_8_120_t ymax = ((m_ymax + m_ymin) * fixed_8_120_t(0.5)) + ysize;
    fixed_8_120_t cr = 0.0, ci = 0.0;
    if (m_SetType == stJulia) {
        cr = m_JuliaCr, ci = m_JuliaCi;
    }

    if (m_zoom > MAX_ZOOM) {
        DrawImageFixedPoint128(pBits, width, height, m_xmin, dx, ymax, -dx, cr, ci);
    }
    else {
        DrawImage(pBits, width, height, m_xmin, dx, ymax, -dx, cr, ci);
    }

    HRESULT hr = image.Save(filename, Gdiplus::ImageFormatPNG);
    if (!SUCCEEDED(hr)) { 
        AfxMessageBox(L"Failed to save image", MB_ICONINFORMATION);
    }
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
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}
