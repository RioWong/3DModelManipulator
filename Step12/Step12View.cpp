// This MFC Samples source code demonstrates using MFC Microsoft Office Fluent User Interface 
// (the "Fluent UI") and is provided only as referential material to supplement the 
// Microsoft Foundation Classes Reference and related electronic documentation 
// included with the MFC C++ library software.  
// License terms to copy, use or distribute the Fluent UI are available separately.  
// To learn more about our Fluent UI licensing program, please visit 
// http://msdn.microsoft.com/officeui.
//
// Copyright (C) Microsoft Corporation
// All rights reserved.

// Step12View.cpp : implementation of the CStep12View class
//

#include "stdafx.h"
// SHARED_HANDLERS can be defined in an ATL project implementing preview, thumbnail
// and search filter handlers and allows sharing of document code with that project.
#ifndef SHARED_HANDLERS
#include "Step12.h"
#endif

#include "Step12Doc.h"
#include "Step12View.h"

#include "MainFrm.h"
#include <iostream>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CStep12View

IMPLEMENT_DYNCREATE(CStep12View, CView)

BEGIN_MESSAGE_MAP(CStep12View, CView)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_KEYDOWN()
END_MESSAGE_MAP()

// CStep12View construction/destruction

CStep12View::CStep12View()
{
	// TODO: add construction code here

}

CStep12View::~CStep12View()
{
}

BOOL CStep12View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: Modify the Window class or styles here by modifying
	//  the CREATESTRUCT cs

	return CView::PreCreateWindow(cs);
}

// CStep12View drawing

void CStep12View::OnDraw(CDC* /*pDC*/)
{
	CStep12Doc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: add draw code for native data here
}

void CStep12View::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CStep12View::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CStep12View diagnostics

#ifdef _DEBUG
void CStep12View::AssertValid() const
{
	CView::AssertValid();
}

void CStep12View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CStep12Doc* CStep12View::GetDocument() const // non-debug version is inline
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CStep12Doc)));
	return (CStep12Doc*)m_pDocument;
}
#endif //_DEBUG


// CStep12View message handlers


int CStep12View::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  在此添加您专用的创建代码
	mOsg = new COsg(m_hWnd);
	//CString tmpMsg("Please open model file...");
	//SendStrToOutputWnd(2, tmpMsg);
	//SendStrToOutputWnd(2, _T("OK"));
	return 0;
}


void CStep12View::OnDestroy()
{
	CView::OnDestroy();

	// TODO: 在此处添加消息处理程序代码
	if(mOsg != 0) {
		delete mOsg;
	}
	WaitForSingleObject(mThreadHandle, 1000);
}


void CStep12View::OnInitialUpdate()
{
	CView::OnInitialUpdate();

	// TODO: 在此添加专用代码和/或调用基类
	CString csFileName = GetDocument()->GetFileName();
    mOsg->initOsg(csFileName.GetString());
    mThreadHandle = (HANDLE)_beginthread(&COsg::Render, 0, mOsg);
	//CString tmpMsg("Open model file successful.");
	//SendStrToOutputWnd(2, tmpMsg);
}


void CStep12View::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值
	mOsg->getViewer()->getEventQueue()->keyPress(nChar);
	
	CMainFrame* pFrame  =  (CMainFrame*)(AfxGetApp()->m_pMainWnd);         
	pFrame->AddStrToOutputWnd(2, "OKKKKK");

	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}
