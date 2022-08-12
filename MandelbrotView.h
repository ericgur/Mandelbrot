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

// MandelbrotView.h : interface of the CMandelbrotView class
//
#pragma once

#include "fixed_point128.h"
using namespace fp128;

typedef fixed_point128<8> fixed_8_120_t;

enum palette_t
{
    palGrey,
    palGradient,
    palVivid,
    palHistogram
};

enum set_type_t
{
    stMandelbrot,
    stJulia,
    stCount
};

enum precision_t
{
    prAuto,
    prDouble,
    prFixedPoint128
};

class CMandelbrotView : public CView
{
protected: // create from serialization only
    CMandelbrotView();
    DECLARE_DYNCREATE(CMandelbrotView)

// Attributes
public:
    CMandelbrotDoc* GetDocument() const;

// Operations
public:

// Overrides
public:
    virtual void OnDraw(CDC* pDC);  // overridden to draw this view
    virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:


// Implementation
public:
    virtual ~CMandelbrotView();
#ifdef _DEBUG
    virtual void AssertValid() const;
    virtual void Dump(CDumpContext& dc) const;
#endif

protected:
    fixed_8_120_t m_xmin, m_xmax, m_ymin, m_ymax;
    double m_zoom;
    int    m_MaxIter;
    BITMAPINFO m_BmpInfo;
    size_t  m_BuffLen;
    COLORREF* m_ColorTable32;
    COLORREF* m_BmpBits;
    int* m_Histogram;
    float m_HsvOffset;
    float* m_Iterations;
    double m_Frequency;
    bool m_NeedToRecompute;
    bool m_SmoothLevel;
    set_type_t m_SetType;
    double m_JuliaCr, m_JuliaCi; // real and imaginary parts of the Julia constant
    palette_t m_PaletteType;
    bool m_IsResizing;
    bool m_AnimatePalette;
    UINT_PTR m_TimerID;
    precision_t m_Precision;

    // Set default values for zoom and coords
    void SetDefaultValues();
    void SetAspectRatio();
    void CreateColorTables();
    void CreateColorTableFromHistogram(float offset);
    void CreateHistogram(const float* pIterations, int width, int height);
    void DrawImageDouble(float* pBits, int width, int height, double x0, double dx, double y0, double dy, double cr = 0, double ci = 0);
    void DrawImageFixedPoint128(float* pBits, int width, int height, const fixed_8_120_t& x0, const fixed_8_120_t& dx, const fixed_8_120_t& y0,
                       const fixed_8_120_t& dy, const fixed_8_120_t& cr, const fixed_8_120_t& ci);
    void CreateDibFromIterations(COLORREF* pBits, const float* pIterations, int width, int height);
    void OnZoomChange(CPoint& point, double zoomMultiplier);

    // Generated message map functions
protected:
    DECLARE_MESSAGE_MAP()
    afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
    afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
    afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
    afx_msg int  OnCreate(LPCREATESTRUCT lpCreateStruct);
    afx_msg void OnSetTypeSelect(UINT nID);
    afx_msg void OnIterationChange(UINT nID);
    afx_msg void OnFileSaveImage(UINT nID);
    afx_msg void OnPaletteChange(UINT nID);
    afx_msg void OnPrecisionSelect(UINT nID);
    afx_msg void OnResetView();
    afx_msg void OnAnimatePalette();
    afx_msg void OnSmoothColorTransitions();
    afx_msg void OnSetTypeChooseJuliaConstant();
    afx_msg void OnEnterSizeMove();
    afx_msg void OnExitSizeMove();
    afx_msg void OnTimer(UINT_PTR nIDEvent);
};

#ifndef _DEBUG  // debug version in MandelbrotView.cpp
inline CMandelbrotDoc* CMandelbrotView::GetDocument() const
   { return reinterpret_cast<CMandelbrotDoc*>(m_pDocument); }
#endif
