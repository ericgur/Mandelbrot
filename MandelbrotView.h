// MandelbrotView.h : interface of the CMandelbrotView class
//

#include "fixed_point128.h"
using namespace fp128;

typedef fixed_point128<8> fixed_8_120_t;

#pragma once
enum set_type_t
{
    stMandelbrot,
    stJulia,
    stCount
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
    float* m_Iterations;
    double m_Frequency;
    bool m_NeedToRedraw;
    bool m_GreyScale;
    bool m_HistogramColoring;
    bool m_SmoothLevel;
    set_type_t m_SetType;
    double m_JuliaCr, m_JuliaCi; // real and imaginary parts of the Julia constant

    // Set default values for zoom and coords
    void SetDefaultValues();
    void SetAspectRatio();
    void CreateColorTables();
    void CreateHsvTable();
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
    afx_msg void OnGreyScale();
    afx_msg void OnFileSaveImage();
    afx_msg void OnHistogramColoring();
    afx_msg void OnResetView();
    afx_msg void OnSetTypeChooseJuliaConstant();
};

#ifndef _DEBUG  // debug version in MandelbrotView.cpp
inline CMandelbrotDoc* CMandelbrotView::GetDocument() const
   { return reinterpret_cast<CMandelbrotDoc*>(m_pDocument); }
#endif
