// MandelbrotView.cpp : implementation of the CMandelbrotView class
//

#include "pch.h"
#include "Mandelbrot.h"
#include "MandelbrotDoc.h"
#include "MandelbrotView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

#define MAX_ZOOM (1ull<<44)
//#define MAX_ZOOM (1ull<<6)


inline COLORREF blendAlpha(COLORREF colora, COLORREF colorb, DWORD alpha)
{
    COLORREF rb1 = (0x100 - alpha) * (colora & 0xFF00FF);
    COLORREF rb2 = alpha * (colorb & 0xFF00FF);
    COLORREF g1 = (0x100 - alpha) * (colora & 0x00FF00);
    COLORREF g2 = alpha * (colorb & 0x00FF00);
    return (((rb1 + rb2) >> 8) & 0xFF00FF) + (((g1 + g2) >> 8) & 0x00FF00);

    //COLORREF rb1 = ((0x100 - alpha) * (colora & 0xFF00FF)) >> 8;
    //COLORREF rb2 = (alpha * (colorb & 0xFF00FF)) >> 8;
    //COLORREF g1 = ((0x100 - alpha) * (colora & 0x00FF00)) >> 8;
    //COLORREF g2 = (alpha * (colorb & 0x00FF00)) >> 8;
    //return ((rb1 + rb2) & 0xFF00FF) + ((g1 + g2) & 0x00FF00);
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
END_MESSAGE_MAP()

// CMandelbrotView construction/destruction

CMandelbrotView::CMandelbrotView()
{
    m_MaxIter = 128;
    m_SmoothLevel = true;

    m_ColorTable32 = NULL;
    m_BmpBits = NULL;
    m_BuffLen = 0;
    
    // set default precision for MPIR
    // TODO: change this depending on zoom level. 64 is good for ~x60
    mpf_set_default_prec(64);

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


void CMandelbrotView::DrawImageMPIR(COLORREF* pBits, int width, int height, const mpf_class& x0, const mpf_class& dx, const mpf_class& y0, const mpf_class& dy)
{
    mpf_class radius = 2.0, radius_sq = radius * radius;

    //create x table
    mpf_class* xTable = new mpf_class[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + dx * i;
    }

#pragma omp parallel for
    for (int l = 0; l < height; ++l) {
        mpf_class y = y0 + (dy * l);
        mpf_class usq = 0, vsq = 0, u = 0, v = 0, x, tmp, uv;
        //point to start of buffer
        COLORREF* pbuff = pBits + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            usq = 0;
            vsq = 0;
            u = 0;
            v = 0;
            x = xTable[k];

            // complex iterative equation is:
            // C(i) = C(i-1) ^ 2 + C(0)
            do {
                // real
                tmp = usq;
                tmp -= vsq;
                tmp += x;

                // imaginary
                //v = 2.0 * (u * v) + y;
                v *= u;
                v *= 2;
                v += y;
                u = tmp;
                vsq = v;
                vsq *= v;
                usq = u;
                usq *= u;

                // check uv vector amplitude is smaller than 2
                tmp = vsq;
                tmp += usq;
                if (tmp >= radius_sq)
                    break;
            } while (++iter < m_MaxIter);

            *(pbuff++) = m_ColorTable32[iter];
        }
    }

    delete[] xTable;
}


void CMandelbrotView::DrawImage(COLORREF* pBits, int width, int height, double x0, double dx, double y0, double dy)
/* Draw the Mandelbrot image on a DIB surface
*     pBits: output DIB surface
*     width: width in pixels
*     height: height in pixels
*     x0: left most coord
*     dx: delta coord between pixels
*     y0: top or bottom most coord. Depending if the image is top down or bottom up
*     dy: delta coord between pixels, negative for top down DIBs
*/
{
    const double radius = 2.0, radius_sq = radius * radius;
    const double LOG2 = log(2.0);

    //create x table
    double* xTable = new double[width];
    for (int i = 0; i < width; ++i) {
        xTable[i] = x0 + (double)(i) * dx;
    }

#pragma omp parallel for
    for (int l = 0; l < height; ++l) {
        double y = y0 + (dy * l);

        //point to start of buffer
        COLORREF* pbuff = pBits + width * l;

        for (int k = 0; k < width; ++k) {
            int iter = 0;
            double usq = 0, vsq = 0, u = 0, v = 0;
            double x = xTable[k];

            // complex iterative equation is:
            // C(i) = C(i-1) ^ 2 + C(0)
            do {
                // real
                double tmp = usq - vsq + x;

                // imaginary
                //v = 2.0 * (u * v) + y;
                v = u * v + u * v + y;
                u = tmp;
                vsq = v * v;
                usq = u * u;
                // check uv vector amplitude is smaller than 2
            } while (vsq + usq < radius_sq && ++iter < m_MaxIter);

            if (m_SmoothLevel && iter < m_MaxIter && iter > 0) {
                // TODO: find working smooth version
                double mu = (double)iter + 1 - (log(log(sqrt(vsq + usq)))) / LOG2;
                DWORD index = (DWORD)floor(mu);
                COLORREF c1 = m_ColorTable32[index];
                COLORREF c2 = m_ColorTable32[index + 1];
                DWORD alpha = (DWORD)(255.0 * (mu - index));
                //COLORREF c1 = m_ColorTable32[iter];
                //COLORREF c2 = m_ColorTable32[iter + 1];
                //DWORD alpha = (DWORD)(128.0 * (sqrt(vsq + usq - radius_sq - x * x))); // 128 = 256 / 2;
                //double fraction = max(0, log(1 + sqrt(vsq + usq) - radius_sq) / LOG2);
                //DWORD alpha = 255 - (DWORD)(255.0 * fraction / radius); 
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
        m_BuffLen = height * width;
        m_BmpBits = (COLORREF*)malloc(m_BuffLen * sizeof(COLORREF));
        m_NeedToRedraw = true;
    }

    mpf_class dx((m_xmax - m_xmin) / width), dy = dx;

    if (m_NeedToRedraw) {
        LARGE_INTEGER time_start, time_end;
        QueryPerformanceCounter(&time_start);

        if (m_zoom > MAX_ZOOM) {
            DrawImageMPIR(m_BmpBits, width, height, m_xmin, dx, m_ymin, dy);
        }
        else {
            DrawImage(m_BmpBits, width, height, m_xmin.get_d(), dx.get_d(), m_ymin.get_d(), dy.get_d());
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


void CMandelbrotView::SetDefaultValues(void)
{
    m_xmax = 2.5;
    m_xmin = -m_xmax;
    m_ymax = m_ymin = 0.0;
    m_zoom = 1;
    SetAspectRatio();
}


void CMandelbrotView::SetAspectRatio(void)
{
    CRect rect;
    GetClientRect(rect);
    //use m_xmin, m_xmax and m_rect to determine m_ymin and m_ymax
    //check if window created
    if (0 == rect.Height())
        return;
    mpf_class ratio = (double)(rect.Height()) / (double)(rect.Width());
    mpf_class ysize((m_xmax - m_xmin) * (ratio / 2.0));
    m_ymin = ((m_ymax + m_ymin) / 2.0) - ysize;
    m_ymax = m_ymin + (2.0 * ysize);
    DebugPrint(L"SetAspectRatio: ysize=%lf, m_ymin=%lf, m_ymax=%lf\n", ysize.get_d(), m_ymin.get_d(), m_ymax.get_d());
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

//zoom in x2
void CMandelbrotView::OnLButtonDown(UINT nFlags, CPoint point)
{
    double zoomMultiplier = 2.0;
    if (nFlags & MK_CONTROL)
        zoomMultiplier = 4.0;

    m_zoom *= zoomMultiplier;
    CRect rect;
    GetClientRect(&rect);

    mpf_class alpha(1.0 - ((double)(point.y) / (double)rect.bottom));
    
    //fix y coords
    mpf_class quarter((m_ymax - m_ymin) / (zoomMultiplier * 2));
    mpf_class center(alpha * m_ymax + (1.0 - alpha) * m_ymin);
    DebugPrint(L"OnLButtonDown: alpha=%lf, quarter=%lf, center=%lf\n", alpha.get_d(), quarter.get_d(), center.get_d());

    m_ymin = center - quarter;
    m_ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (m_xmax - m_xmin) / (zoomMultiplier * 2.0);
    center = alpha * m_xmax + (1.0 - alpha) * m_xmin;
    m_xmin = center - quarter;
    m_xmax = center + quarter;

    m_NeedToRedraw = true;
    Invalidate(FALSE);

    // CView::OnLButtonDown(nFlags, point);
}


//zoom out x2
void CMandelbrotView::OnRButtonDown(UINT nFlags, CPoint point)
{
    double zoomMultiplier = 2.0;
    if (nFlags & MK_CONTROL)
        zoomMultiplier = 4.0;

    CRect rect;
    GetClientRect(&rect);

    m_zoom /= zoomMultiplier;

    mpf_class alpha(1.0 - ((double)(point.y) / (double)rect.bottom));
    //fix y coords
    mpf_class quarter((m_ymax - m_ymin) * zoomMultiplier / 2.0);
    mpf_class center(alpha * m_ymax + (1.0 - alpha) * m_ymin);
    m_ymin = center - quarter;
    m_ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (m_xmax - m_xmin) * zoomMultiplier / 2.0;
    center = alpha * m_xmax + (1.0 - alpha) * m_xmin;
    m_xmin = center - quarter;
    m_xmax = center + quarter;

    m_NeedToRedraw = true;
    Invalidate(FALSE);

    //CView::OnRButtonDown(nFlags, point);
}


//reset to default
void CMandelbrotView::OnMButtonDown(UINT nFlags, CPoint point)
{
    SetDefaultValues();
    m_NeedToRedraw = true;
    Invalidate(FALSE);

    CView::OnMButtonDown(nFlags, point);
}


void CMandelbrotView::CreateColorTables(void)
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
    //for (UINT i = ID_ITERATIONS; i <= ID_ITERATIONS_LAST; ++i) {
    //    menu->CheckMenuItem(i, (i == nID) ? MF_CHECKED : MF_UNCHECKED);
    //}
    CreateColorTables();
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}


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

    mpf_class dx((m_xmax - m_xmin) / width);

    if (m_zoom > MAX_ZOOM) {
        DrawImageMPIR(pBits, width, height, m_xmin, dx, m_ymin, -dx);
        return;
    }
    else {
        DrawImage(pBits, width, height, m_xmin.get_d(), dx.get_d(), m_ymin.get_d(), -dx.get_d());
    }

    HRESULT hr = image.Save(filename, Gdiplus::ImageFormatPNG);
    if (!SUCCEEDED(hr)) { 
        AfxMessageBox(L"Failed to save image", MB_ICONINFORMATION);
    }
}
