
// Step11View.h : CStep11View ��Ľӿ�
//

#pragma once
#include "Osg.h"


class CStep11View : public CView
{
protected: // �������л�����
	CStep11View();
	DECLARE_DYNCREATE(CStep11View)

// ����
public:
	CStep11Doc* GetDocument() const;

// ����
public:

// ��д
public:
	virtual void OnDraw(CDC* pDC);  // ��д�Ի��Ƹ���ͼ
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:

// ʵ��
public:
	virtual ~CStep11View();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	COsg * mOsg;
	HANDLE mThreadHandle;

// ���ɵ���Ϣӳ�亯��
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
public:
	afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);
	afx_msg void OnDestroy();
	virtual void OnInitialUpdate();
	afx_msg void OnKeyDown(UINT nChar, UINT nRepCnt, UINT nFlags);
	afx_msg BOOL OnEraseBkgnd(CDC* pDC);
};

#ifndef _DEBUG  // Step11View.cpp �еĵ��԰汾
inline CStep11Doc* CStep11View::GetDocument() const
   { return reinterpret_cast<CStep11Doc*>(m_pDocument); }
#endif

