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

// MandelbrotDoc.cpp : implementation of the CMandelbrotDoc class
//

#include "pch.h"
#include "Mandelbrot.h"
#include "MandelbrotDoc.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CMandelbrotDoc

IMPLEMENT_DYNCREATE(CMandelbrotDoc, CDocument)

BEGIN_MESSAGE_MAP(CMandelbrotDoc, CDocument)
END_MESSAGE_MAP()


// CMandelbrotDoc construction/destruction

CMandelbrotDoc::CMandelbrotDoc()
{
}


CMandelbrotDoc::~CMandelbrotDoc()
{
}


BOOL CMandelbrotDoc::OnNewDocument()
{
    if (!CDocument::OnNewDocument())
        return FALSE;

    // add reinitialization code here
    // (SDI documents will reuse this document)

    return TRUE;
}


// CMandelbrotDoc serialization
void CMandelbrotDoc::Serialize(CArchive& ar)
{
    if (ar.IsStoring()) {
        // TODO: add storing code here
    }
    else {
        // TODO: add loading code here
    }
}


// CMandelbrotDoc diagnostics
#ifdef _DEBUG
void CMandelbrotDoc::AssertValid() const
{
    CDocument::AssertValid();
}

void CMandelbrotDoc::Dump(CDumpContext& dc) const
{
    CDocument::Dump(dc);
}
#endif //_DEBUG


// CMandelbrotDoc commands
