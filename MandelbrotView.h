// MandelbrotView.h : interface of the CMandelbrotView class
//


#pragma once


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
    double m_xmin, m_xmax, m_ymin, m_ymax;
    double m_zoom;
    size_t m_MaxIter;
    BITMAPINFO m_BmpInfo;
    int m_BuffLen;
    COLORREF* m_ColorTable32;
    COLORREF* m_BmpBits;
    double m_Frequency;
    bool m_NeedToRedraw;
    bool m_GreyScale;
    bool m_SmoothLevel;

    // Set default values for zoom and coords
    void SetDefaultValues(void);
    void SetAspectRatio(void);
    void CreateColorTables(void);
    void DrawImage(COLORREF* pBits, int width, int height, double x0, double dx, double y0, double dy);
// Generated message map functions
protected:
    DECLARE_MESSAGE_MAP()
    afx_msg void OnLButtonDown(UINT nFlags, CPoint point);
    afx_msg void OnRButtonDown(UINT nFlags, CPoint point);
    afx_msg void OnMButtonDown(UINT nFlags, CPoint point);
    afx_msg int  OnCreate(LPCREATESTRUCT lpCreateStruct);
    afx_msg void OnIterationChange(UINT nID);
    afx_msg void OnGreyScale();
    afx_msg void OnFileSaveImage();
};

#ifndef _DEBUG  // debug version in MandelbrotView.cpp
inline CMandelbrotDoc* CMandelbrotView::GetDocument() const
   { return reinterpret_cast<CMandelbrotDoc*>(m_pDocument); }
#endif

