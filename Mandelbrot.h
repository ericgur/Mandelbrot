// Mandelbrot.h : main header file for the Mandelbrot application
//
#pragma once

#ifndef __AFXWIN_H__
    #error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"       // main symbols


// CMandelbrotApp:
// See Mandelbrot.cpp for the implementation of this class
//

class CMandelbrotApp : public CWinApp
{
public:
    CMandelbrotApp();


// Overrides
public:
    virtual BOOL InitInstance();

// Implementation
    afx_msg void OnAppAbout();
    DECLARE_MESSAGE_MAP()
};

extern CMandelbrotApp theApp;