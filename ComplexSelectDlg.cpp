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

	return TRUE;  // return TRUE unless you set the focus to a control  
}

void ComplexSelectDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//DDX_Control(pDX, IDC_EDIT_IMAG, m_Imag);
	//DDX_Control(pDX, IDC_EDIT_REAL, m_Real);
	DDX_Text(pDX, IDC_EDIT_REAL, m_Real);
	DDX_Text(pDX, IDC_EDIT_IMAG, m_Imag);
}


BEGIN_MESSAGE_MAP(ComplexSelectDlg, CDialog)
END_MESSAGE_MAP()


// ComplexSelectDlg message handlers
