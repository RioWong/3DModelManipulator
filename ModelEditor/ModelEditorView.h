
// ModelEditorView.h : interface of the CModelEditorView class
//

#pragma once


class CModelEditorView : public CView
{
protected: // create from serialization only
	CModelEditorView();
	DECLARE_DYNCREATE(CModelEditorView)

// Attributes
public:
	CModelEditorDoc* GetDocument() const;

// Operations
public:

// Overrides
public:
	virtual void OnDraw(CDC* pDC);  // overridden to draw this view
	virtual BOOL PreCreateWindow(CREATESTRUCT& cs);
protected:

// Implementation
public:
	virtual ~CModelEditorView();
#ifdef _DEBUG
	virtual void AssertValid() const;
	virtual void Dump(CDumpContext& dc) const;
#endif

protected:

// Generated message map functions
protected:
	afx_msg void OnFilePrintPreview();
	afx_msg void OnRButtonUp(UINT nFlags, CPoint point);
	afx_msg void OnContextMenu(CWnd* pWnd, CPoint point);
	DECLARE_MESSAGE_MAP()
};

#ifndef _DEBUG  // debug version in ModelEditorView.cpp
inline CModelEditorDoc* CModelEditorView::GetDocument() const
   { return reinterpret_cast<CModelEditorDoc*>(m_pDocument); }
#endif

