#pragma once
#include <osgViewer/Viewer>
#include <osg/Node>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <iostream>
#include <osgManipulator/TranslateAxisDragger>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/NodeVisitor>
#include <osg/Point>

#include <osgGA/TrackballManipulator>
#include <osgManipulator/TrackballDragger>
#include <osgGA/KeySwitchMatrixManipulator>
#include <osgViewer/api/win32/GraphicsWindowWin32>

#include <iostream>

class COsg
{
public:
	COsg(HWND hWnd);
	~COsg(void);
	void initOsg(std::string fileName);
	osgViewer::Viewer * getViewer(void);
	static void Render(void* ptr);

private:
	void initCameraConfig(void);
	void initManipulators(void);

private:
	HWND m_hWnd;
	std::string mFileName;
	osgViewer::Viewer * mViewer;
	osg::ref_ptr<osg::Group> mRoot;
	osg::ref_ptr<osg::Node> mModel;
	osg::ref_ptr<osgGA::TrackballManipulator> trackball;
	osg::ref_ptr<osgGA::KeySwitchMatrixManipulator> keyswitchManipulator;
};

