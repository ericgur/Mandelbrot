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

// ComplexSelectDlg.cpp : implementation file
//

#include "pch.h"
#include "Mandelbrot.h"
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
    presets.AddString(L"0.0, -0.8");
    presets.AddString(L"0.45, 0.1428");
    presets.AddString(L"-0.70176, -0.3842");
    presets.AddString(L"-0.75, 0.11");
    presets.AddString(L"-0.1, 0.651");
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
    auto index = presets.GetCurSel();
    CString txt, s;
    presets.GetText(index, txt);
    int pos = 0;
    real = txt.Tokenize(L", ", pos);
    imag = txt.Tokenize(L", ", pos);
    UpdateData(FALSE);
}
