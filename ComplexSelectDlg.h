#pragma once
#include "afxdialogex.h"

// ComplexSelectDlg dialog

class ComplexSelectDlg : public CDialog
{
    DECLARE_DYNAMIC(ComplexSelectDlg)

public:
    ComplexSelectDlg(CWnd* pParent = nullptr);   // standard constructor
    virtual ~ComplexSelectDlg();
    virtual BOOL OnInitDialog();
// Dialog Data
#ifdef AFX_DESIGN_TIME
    enum { IDD = IDD_COMPLEX_SELECT };
#endif

protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

    DECLARE_MESSAGE_MAP()
public:
    CString real;
    CString imag;
    CListBox presets;
    afx_msg void OnLbnSelchangeListPresets();
};
