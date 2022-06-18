// MandelbrotView.cpp : implementation of the CMandelbrotView class
//

#include "pch.h"
#include "Mandelbrot.h"

#include "MandelbrotDoc.h"
#include "MandelbrotView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

static const int s_iterCountTable[] =
{
    64,
    96,
    128,
    192,
    256,
    384,
    512
};

// CMandelbrotView

IMPLEMENT_DYNCREATE(CMandelbrotView, CView)

BEGIN_MESSAGE_MAP(CMandelbrotView, CView)
ON_WM_LBUTTONDOWN()
ON_WM_RBUTTONDOWN()
ON_WM_MBUTTONDOWN()
ON_WM_CREATE()
ON_COMMAND(ID_VIEW_GREYSCALE, OnGreyScale)
ON_COMMAND_RANGE(ID_ITERATIONS, ID_ITERATIONS_LAST, OnIterationChange)
END_MESSAGE_MAP()

// CMandelbrotView construction/destruction

CMandelbrotView::CMandelbrotView() 
{
    m_MaxIter = 128;
    m_ColorTable32 = NULL;
    m_BmpBits = NULL;
    m_BuffLen = 0;

    //fill bitmap header
    memset(&m_BmpInfo, 0, sizeof(m_BmpInfo));
    m_BmpInfo.bmiHeader.biSize          = sizeof(BITMAPINFOHEADER);
    m_BmpInfo.bmiHeader.biPlanes        = 1; 
    m_BmpInfo.bmiHeader.biBitCount      = 32;
    m_BmpInfo.bmiHeader.biCompression   = BI_RGB;
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


void CMandelbrotView::OnDraw(CDC* pDC)
{
    SetAspectRatio();

    wchar_t message[256];

    CRect rect;
    GetClientRect(rect);
    const int width = rect.Width(), height = rect.Height();

    //deallocate on size change
    if(m_BmpInfo.bmiHeader.biHeight != height || m_BmpInfo.bmiHeader.biWidth != width)
    {
        free(m_BmpBits);
        m_BmpBits = NULL;
    }

    //allocate new bitmap if needed
    if (NULL == m_BmpBits)
    {
        m_BmpInfo.bmiHeader.biHeight = height;
        m_BmpInfo.bmiHeader.biWidth  = width;
        m_BuffLen = height * width;
        m_BmpBits = (COLORREF*)malloc(m_BuffLen * sizeof(COLORREF));
        m_NeedToRedraw = true;
    }
    
    const double dx = (xmax-xmin) / width, dy = (ymax-ymin) / height;

    if (xmin + dx == xmin || ymin + dy == ymin)
    {
        AfxMessageBox(L"Maximum precision reached :)\nPlease zoom out", MB_ICONINFORMATION);
        m_zoom *= 4.0;
        m_NeedToRedraw = false;
    }

    if(!m_NeedToRedraw)
        return;

    LARGE_INTEGER time_start, time_end;
    QueryPerformanceCounter(&time_start);
    

    //create x table
    double* xTable = new double[width];
    xTable[0] = xmin;
    for(int i=1; i<width; ++i)
    {
        xTable[i] = xmin + ((double)(i) * dx);
    }

#pragma omp parallel for

    for (int l = 0; l < height; ++l)
    {
        double y = ymin + (dy * l);

        //point to start of buffer
        COLORREF* pbuff = m_BmpBits + width * l;

        for (int k = 0; k < width; ++k)
        {
            int color = 0;
            double usq = 0, vsq = 0, u = 0, v = 0;
            double x = xTable[k];

            // complex iterative equation is:
            // C(i) = C(i-1)^2 + C(0)
            do
            {
                // real
                double tmp = usq - vsq + x;

                // imaginary
                //v = 2.0*(u*v)+ y;
                v = u*v + u*v + y;
                u = tmp;
                vsq = v*v;
                usq = u*u;
                // check uv vector amplitude is smaller than 2
            } while (vsq+usq < 4.0 && ++color < m_MaxIter);

            *(pbuff++) = m_ColorTable32[color];
        }        
    }

    //all done
    QueryPerformanceCounter(&time_end);

    DWORD totalTime = DWORD(1000.0 * (time_end.QuadPart - time_start.QuadPart) / m_Frequency);
    delete[] xTable;

    SetDIBitsToDevice((HDC)(*pDC), 0, 0, width, height, 0, 0, 0, height, m_BmpBits, &m_BmpInfo, DIB_RGB_COLORS);

    swprintf(message, _countof(message), (m_zoom >= 1.0) ? L"Zoom x%0.0lf (%ims)" : L"Zoom x%0.5f (%ims)", m_zoom, totalTime);

    ((CFrameWnd*)AfxGetMainWnd())->SetWindowText( message );
    m_NeedToRedraw = false;
}


void CMandelbrotView::SetDefaultValues(void)
{
    xmax = 2.5f;
    xmin = -xmax;
    m_zoom = 1;
    SetAspectRatio();
}


void CMandelbrotView::SetAspectRatio(void)
{
    CRect rect;
    GetClientRect(rect);
    //use xmin, xmax and m_rect to determine ymin and ymax
    //check if window created
    if (0 == rect.Height())
        return;
    double ratio = (double)(rect.Height()) / (double)(rect.Width());
    double ysize = (xmax - xmin) * (ratio / 2.0);
    ymin = ((ymax + ymin) / 2.0) - ysize;
    ymax = ymin + (2.0 * ysize);
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
    double quarter, center, alpha;
    CRect rect;
    GetClientRect( &rect );

    alpha = 1.0 - ((double)(point.y) / (double)rect.bottom);
    //fix y coords
    quarter = (ymax - ymin) / 4.0;
    center = alpha * ymax + (1.0 - alpha) * ymin;
    ymin = center - quarter;
    ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (xmax - xmin) / 4.0;
    center = alpha * xmax + (1.0 - alpha) * xmin;
    xmin = center - quarter;
    xmax = center + quarter;

    m_zoom *= 2.0;
    m_NeedToRedraw = true;
    Invalidate(FALSE);

    CView::OnLButtonDown(nFlags, point);
}


//zoom out x2
void CMandelbrotView::OnRButtonDown(UINT nFlags, CPoint point)
{
    double quarter, center, alpha;
    CRect rect;
    GetClientRect( &rect );

    alpha = 1.0 - ((double)(point.y) / (double)rect.bottom);
    //fix y coords
    quarter = (ymax - ymin);
    center = alpha * ymax + (1.0 - alpha) * ymin;
    ymin = center - quarter;
    ymax = center + quarter;

    //fix x coords
    alpha = (double)(point.x) / (double)rect.right;
    quarter = (xmax - xmin);
    center = alpha * xmax + (1.0 - alpha) * xmin;
    xmin = center - quarter;
    xmax = center + quarter;

    m_zoom /= 2.0;

    m_NeedToRedraw = true;
    Invalidate(FALSE);

    CView::OnRButtonDown(nFlags, point);
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

    m_ColorTable32 = new COLORREF[m_MaxIter+2];
    if (m_GreyScale)
    {
        for (size_t i = 1; i <= m_MaxIter; ++i)
        {
            int c = 255 - (int)(255.0f * (float)i / (float)m_MaxIter);
            m_ColorTable32[i] = RGB(c, c, c);
        }
    }
    else
    {
        for (size_t i = 1; i <= m_MaxIter; ++i)
        {
            size_t c = m_MaxIter - i;
            m_ColorTable32[i] = RGB(c * 12 & 255, c * 16 & 255, c * 5 & 255);
        }
    }

    m_ColorTable32[0] = RGB(255, 255, 255);
}


int CMandelbrotView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
    if (CView::OnCreate(lpCreateStruct) == -1)
        return -1;

    SetDefaultValues();

    return 0;
}


void CMandelbrotView::OnIterationChange(UINT nID)
{ 
    m_MaxIter = s_iterCountTable[nID - ID_ITERATIONS];

    CMenu *menu = AfxGetMainWnd()->GetMenu();

    for (int i = ID_ITERATIONS; i <= ID_ITERATIONS + 6; ++i)
    {
        if((UINT)i == nID)
            menu->CheckMenuItem(i, MF_CHECKED);
        else
            menu->CheckMenuItem(i, MF_UNCHECKED);
    }
    CreateColorTables();
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}


void CMandelbrotView::OnGreyScale()
{
    CMenu* menu = AfxGetMainWnd()->GetMenu();

    int state = menu->GetMenuState(ID_VIEW_GREYSCALE, MF_BYCOMMAND);

    if (state & MF_CHECKED)
    {
        menu->CheckMenuItem(ID_VIEW_GREYSCALE, MF_UNCHECKED | MF_BYCOMMAND);
        m_GreyScale = false;
    }
    else
    {
        menu->CheckMenuItem(ID_VIEW_GREYSCALE, MF_CHECKED | MF_BYCOMMAND);
        m_GreyScale = true;
    }
    
    CreateColorTables();
    m_NeedToRedraw = true;
    Invalidate(FALSE);
}