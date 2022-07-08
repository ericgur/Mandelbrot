// ComplexSelectDlg.cpp : implementation file
//

#include "pch.h"
#include "Mandelbrot.h"
#include "afxdialogex.h"
#include "ComplexSelectDlg.h"


// ComplexSelectDlg dialog

IMPLEMENT_DYNAMIC(ComplexSelectDlg, CDialog)

ComplexSelectDlg::ComplexSelectDlg(CWnd* pParent /*=nullptr*/)
    : CDialog(IDD_COMPLEX_SELECT, pParent)
{

}

ComplexSelectDlg::~ComplexSelectDlg()
{
}


BOOL ComplexSelectDlg::OnInitDialog()
{
    CDialog::OnInitDialog();
    presets.AddString(L"0.285, 0.01");
    presets.AddString(L"-0.7269, 0.1889");
    presets.AddString(L"-0.8, 0.156");
    presets.AddString(L"-0.4, 0.6");

    return TRUE;  // return TRUE unless you set the focus to a control  
}

void ComplexSelectDlg::DoDataExchange(CDataExchange* pDX)
{
    CDialog::DoDataExchange(pDX);
    DDX_Text(pDX, IDC_EDIT_REAL, real);
    DDX_Text(pDX, IDC_EDIT_IMAG, imag);
    DDX_Control(pDX, IDC_LIST_PRESETS, presets);
}


BEGIN_MESSAGE_MAP(ComplexSelectDlg, CDialog)
    ON_LBN_SELCHANGE(IDC_LIST_PRESETS, &ComplexSelectDlg::OnLbnSelchangeListPresets)
END_MESSAGE_MAP()


// ComplexSelectDlg message handlers


void ComplexSelectDlg::OnLbnSelchangeListPresets()
{
    // TODO: Add your control notification handler code here
    auto index = presets.GetCurSel();
    CString txt, s;
    presets.GetText(index, txt);
    int pos = 0;
    //= txt.Tokenize(L", ", pos);
    real = txt.Tokenize(L", ", pos);
    imag = txt.Tokenize(L", ", pos);
    UpdateData(FALSE);
}
