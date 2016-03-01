#pragma once
#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osgUtil/Optimizer>
#include <osg/Node>
#include <osgManipulator/TrackballDragger>
#include <osgManipulator/TranslateAxisDragger>
class COSG
{
public:
	COSG(HWND hWnd);
	~COSG(void);

	void initOSG(char * fileName);
	void InitManipulators(void);
	void InitSceneGraph(void);
	void InitCameraConfig(void);
	void SetupWindow(void);
	void SetupCamera(void);
	void PreFrameUpdate(void);
	void PostFrameUpdate(void);
	void Done(bool value) { mDone = value; }
	bool Done(void) { return mDone; }
	static void Render(void* ptr);

	osgViewer::Viewer* getViewer() { return mViewer; }

private:
	bool mDone;
	std::string m_ModelName;
	HWND mHWnd;
	osgViewer::Viewer* mViewer;
	osg::ref_ptr<osg::Group> mRoot;
	osg::ref_ptr<osg::Node> mModel;

};

