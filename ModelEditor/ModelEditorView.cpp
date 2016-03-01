
// ModelEditorView.cpp : implementation of the CModelEditorView class
//

#include "stdafx.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "ModelEditor.h"
#endif

#include "ModelEditorDoc.h"
#include "ModelEditorView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CModelEditorView

IMPLEMENT_DYNCREATE(CModelEditorView, CView)

BEGIN_MESSAGE_MAP(CModelEditorView, CView)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
END_MESSAGE_MAP()

// CModelEditorView construction/destruction

CModelEditorView::CModelEditorView()
{
	// TODO: add construction code here

}

CModelEditorView::~CModelEditorView()
{
}

BOOL CModelEditorView::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CModelEditorView drawing

void CModelEditorView::OnDraw(CDC* /*pDC*/)
{
	CModelEditorDoc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: add draw code for native data here
}

void CModelEditorView::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CModelEditorView::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CModelEditorView diagnostics

#ifdef _DEBUG
void CModelEditorView::AssertValid() const
{
	CView::AssertValid();
}

void CModelEditorView::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CModelEditorDoc* CModelEditorView::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CModelEditorDoc)));
	return (CModelEditorDoc*)m_pDocument;
}
#endif //_DEBUG


// CModelEditorView message handlers
