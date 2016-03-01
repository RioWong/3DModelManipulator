#include "StdAfx.h"
#include "OSG.h"


COSG::COSG(HWND hWnd) : mHWnd(hWnd)
{
}


COSG::~COSG(void)
{
	mViewer->setDone(true);
	Sleep(1000);
	mViewer->stopThreading();
	delete mViewer;
}

void initOSG(char * fileName)
{
	
}
