
// Step11View.cpp : CStep11View 类的实现
//

#include "stdafx.h"
// SHARED_HANDLERS 可以在实现预览、缩略图和搜索筛选器句柄的
// ATL 项目中进行定义，并允许与该项目共享文档代码。
#ifndef SHARED_HANDLERS
#include "Step11.h"
#endif

#include "Step11Doc.h"
#include "Step11View.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CStep11View

IMPLEMENT_DYNCREATE(CStep11View, CView)

BEGIN_MESSAGE_MAP(CStep11View, CView)
	ON_WM_CONTEXTMENU()
	ON_WM_RBUTTONUP()
	ON_WM_CREATE()
	ON_WM_DESTROY()
	ON_WM_KEYDOWN()
	ON_WM_ERASEBKGND()
END_MESSAGE_MAP()

// CStep11View 构造/析构

CStep11View::CStep11View()
{
	// TODO: 在此处添加构造代码

}

CStep11View::~CStep11View()
{
}

BOOL CStep11View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: 在此处通过修改
	//  CREATESTRUCT cs 来修改窗口类或样式

	return CView::PreCreateWindow(cs);
}

// CStep11View 绘制

void CStep11View::OnDraw(CDC* /*pDC*/)
{
	CStep11Doc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: 在此处为本机数据添加绘制代码
}

void CStep11View::OnRButtonUp(UINT /* nFlags */, CPoint point)
{
	ClientToScreen(&point);
	OnContextMenu(this, point);
}

void CStep11View::OnContextMenu(CWnd* /* pWnd */, CPoint point)
{
#ifndef SHARED_HANDLERS
	theApp.GetContextMenuManager()->ShowPopupMenu(IDR_POPUP_EDIT, point.x, point.y, this, TRUE);
#endif
}


// CStep11View 诊断

#ifdef _DEBUG
void CStep11View::AssertValid() const
{
	CView::AssertValid();
}

void CStep11View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CStep11Doc* CStep11View::GetDocument() const // 非调试版本是内联的
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CStep11Doc)));
	return (CStep11Doc*)m_pDocument;
}
#endif //_DEBUG


// CStep11View 消息处理程序


int CStep11View::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  在此添加您专用的创建代码
	mOsg = new COsg(m_hWnd);
	return 1;
}


void CStep11View::OnDestroy()
{
	CView::OnDestroy();

	// TODO: 在此处添加消息处理程序代码
	if(mOsg != 0) {
		delete mOsg;
	}
	WaitForSingleObject(mThreadHandle, 1000);
}


void CStep11View::OnInitialUpdate()
{
	CView::OnInitialUpdate();

	// TODO: 在此添加专用代码和/或调用基类
	CString csFileName = GetDocument()->GetFileName();
	mOsg->initOsg(csFileName.GetString());
	//mOsg->initOsg("glider.osg");
	mThreadHandle = (HANDLE)_beginthread(&COsg::Render, 0, mOsg);
}


void CStep11View::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值
	mOsg->getViewer()->getEventQueue()->keyPress(nChar);

	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}


BOOL CStep11View::OnEraseBkgnd(CDC* pDC)
{
	// TODO: 在此添加消息处理程序代码和/或调用默认值

	//return CView::OnEraseBkgnd(pDC);
	/* Do nothing, to avoid flashing on MSW */
	return true;
}
