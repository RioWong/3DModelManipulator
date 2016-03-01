
// Step11View.cpp : CStep11View ���ʵ��
//

#include "stdafx.h"
// SHARED_HANDLERS ������ʵ��Ԥ��������ͼ������ɸѡ�������
// ATL ��Ŀ�н��ж��壬�����������Ŀ�����ĵ����롣
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

// CStep11View ����/����

CStep11View::CStep11View()
{
	// TODO: �ڴ˴���ӹ������

}

CStep11View::~CStep11View()
{
}

BOOL CStep11View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: �ڴ˴�ͨ���޸�
	//  CREATESTRUCT cs ���޸Ĵ��������ʽ

	return CView::PreCreateWindow(cs);
}

// CStep11View ����

void CStep11View::OnDraw(CDC* /*pDC*/)
{
	CStep11Doc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	// TODO: �ڴ˴�Ϊ����������ӻ��ƴ���
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


// CStep11View ���

#ifdef _DEBUG
void CStep11View::AssertValid() const
{
	CView::AssertValid();
}

void CStep11View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CStep11Doc* CStep11View::GetDocument() const // �ǵ��԰汾��������
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CStep11Doc)));
	return (CStep11Doc*)m_pDocument;
}
#endif //_DEBUG


// CStep11View ��Ϣ�������


int CStep11View::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
	if (CView::OnCreate(lpCreateStruct) == -1)
		return -1;

	// TODO:  �ڴ������ר�õĴ�������
	mOsg = new COsg(m_hWnd);
	return 1;
}


void CStep11View::OnDestroy()
{
	CView::OnDestroy();

	// TODO: �ڴ˴������Ϣ����������
	if(mOsg != 0) {
		delete mOsg;
	}
	WaitForSingleObject(mThreadHandle, 1000);
}


void CStep11View::OnInitialUpdate()
{
	CView::OnInitialUpdate();

	// TODO: �ڴ����ר�ô����/����û���
	CString csFileName = GetDocument()->GetFileName();
	mOsg->initOsg(csFileName.GetString());
	//mOsg->initOsg("glider.osg");
	mThreadHandle = (HANDLE)_beginthread(&COsg::Render, 0, mOsg);
}


void CStep11View::OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags)
{
	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ
	mOsg->getViewer()->getEventQueue()->keyPress(nChar);

	CView::OnKeyDown(nChar, nRepCnt, nFlags);
}


BOOL CStep11View::OnEraseBkgnd(CDC* pDC)
{
	// TODO: �ڴ������Ϣ�����������/�����Ĭ��ֵ

	//return CView::OnEraseBkgnd(pDC);
	/* Do nothing, to avoid flashing on MSW */
	return true;
}
