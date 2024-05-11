/***********************************************************************************
    MIT License

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

// MainFrm.cpp : implementation of the CMainFrame class
//

#include "pch.h"
#include "framework.h"
#include "Mandelbrot.h"
#include "MainFrm.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMainFrame

IMPLEMENT_DYNCREATE(CMainFrame, CFrameWnd)

BEGIN_MESSAGE_MAP(CMainFrame, CFrameWnd)
    ON_WM_CREATE()
    ON_WM_ENTERSIZEMOVE()
    ON_WM_EXITSIZEMOVE()
END_MESSAGE_MAP()

static UINT indicators[] =
{
    ID_SEPARATOR,           // status line indicator
    ID_INDICATOR_CAPS,
    ID_INDICATOR_NUM,
    ID_INDICATOR_SCRL,
};


// CMainFrame construction/destruction

CMainFrame::CMainFrame()
{
    // add member initialization code here
}

CMainFrame::~CMainFrame()
{
}


int CMainFrame::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
    if (CFrameWnd::OnCreate(lpCreateStruct) == -1)
        return -1;

    //if (!m_wndToolBar.CreateEx(this, TBSTYLE_FLAT, WS_CHILD | WS_VISIBLE | CBRS_TOP
    //    | CBRS_GRIPPER | CBRS_TOOLTIPS | CBRS_FLYBY | CBRS_SIZE_DYNAMIC) ||
    //    !m_wndToolBar.LoadToolBar(IDR_MAINFRAME))
    //{
    //    TRACE0("Failed to create toolbar\n");
    //    return -1;      // fail to create
    //}

    if (!m_wndStatusBar.Create(this) ||
        !m_wndStatusBar.SetIndicators(indicators, sizeof(indicators) / sizeof(UINT))) {
        TRACE0("Failed to create status bar\n");
        return -1;      // fail to create
    }

    //hide the pretty much unused status bar.
    m_wndStatusBar.ShowWindow(SW_HIDE);

    // Delete these three lines if you don't want the toolbar to be dockable
    //m_wndToolBar.EnableDocking(CBRS_ALIGN_ANY);
    //EnableDocking(CBRS_ALIGN_ANY);
    //DockControlBar(&m_wndToolBar);

    CMenu* menu = GetMenu();
    menu->CheckMenuRadioItem(ID_ITERATIONS, ID_ITERATIONS_LAST, ID_ITERATIONS + 2, MF_BYCOMMAND);
    menu->CheckMenuRadioItem(ID_SETTYPE_MANDELBROT, ID_SETTYPE_JULIA, ID_SETTYPE_MANDELBROT, MF_BYCOMMAND);
    menu->CheckMenuRadioItem(ID_VIEW_GREYSCALE, ID_VIEW_HISTOGRAMCOLORING, ID_VIEW_GRADIENTS, MF_BYCOMMAND);
    menu->CheckMenuItem(ID_VIEW_SMOOTHCOLORTRANSITION, MF_CHECKED | MF_BYCOMMAND);
    menu->CheckMenuRadioItem(ID_PRECISION_AUTO, ID_PRECISION_FIXEDPOINT128, ID_PRECISION_AUTO, MF_BYCOMMAND);
    return 0;
}

BOOL CMainFrame::PreCreateWindow(CREATESTRUCT& cs)
{
    if (!CFrameWnd::PreCreateWindow(cs))
        return FALSE;
    // Modify the Window class or styles here by modifying
    //  the CREATESTRUCT cs

    return TRUE;
}


// CMainFrame diagnostics

#ifdef _DEBUG
void CMainFrame::AssertValid() const
{
    CFrameWnd::AssertValid();
}

void CMainFrame::Dump(CDumpContext& dc) const
{
    CFrameWnd::Dump(dc);
}

#endif //_DEBUG


// CMainFrame message handlers

