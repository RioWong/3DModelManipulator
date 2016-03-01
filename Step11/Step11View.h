
// Step11View.h : CStep11View 类的接口
//

#pragma once
#include "Osg.h"


class CStep11View : public CView
{
protected: // 仅从序列化创建
	CStep11View();
	DECLARE_DYNCREATE(CStep11View)

// 特性
public:
	CStep11Doc* GetDocument() const;

// 操作
public:

// 重写
public:
	virtual void OnDraw(CDC* pDC);  // 重写以绘制该视图
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:

// 实现
public:
	virtual ~CStep11View();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:
	COsg * mOsg;
	HANDLE mThreadHandle;

// 生成的消息映射函数
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

#ifndef _DEBUG  // Step11View.cpp 中的调试版本
inline CStep11Doc* CStep11View::GetDocument() const
   { return reinterpret_cast<CStep11Doc*>(m_pDocument); }
#endif

