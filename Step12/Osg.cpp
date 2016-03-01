#include "StdAfx.h"
#include "Osg.h"


COsg::COsg(HWND hWnd) : m_hWnd(hWnd)
{
}


COsg::~COsg(void)
{
	mViewer->setDone(true);
	Sleep(1000);
	mViewer->stopThreading();
	delete mViewer;
}

void COsg::initOsg(std::string fileName)
{
	/*
	if(fileName == "") {
		mFileName = "glider.osg";
	} else {
		//ɾ��mRoot->getChild(0);
		mFileName = fileName;
	}*/
	mFileName = fileName;
	//mFileName = "glider.osg";
	initManipulators();
	mRoot = new osg::Group();
	//mModel = osgDB::readNodeFile(mFileName);
	osg::ref_ptr<osg::Node> ply2Model = readPly2(fileName);
	osg::ref_ptr<osg::Node> rrModel = readPwn("d:\\osg2014\\model\\bdf2.rpn");

	mRoot->addChild(ply2Model.get());
	mRoot->addChild(rrModel.get());
	
	initCameraConfig();
	osgDB::writeNodeFile(*mRoot, "d:\\osg2014\\step12.osg");
}
osgViewer::Viewer * COsg::getViewer(void)
{
	return mViewer;
}
void COsg::initManipulators(void)
{
    // Create a trackball manipulator
    trackball = new osgGA::TrackballManipulator();

    // Create a Manipulator Switcher
    keyswitchManipulator = new osgGA::KeySwitchMatrixManipulator;

    // Add our trackball manipulator to the switcher
    keyswitchManipulator->addMatrixManipulator( '1', "Trackball", trackball.get());

    // Init the switcher to the first manipulator (in this case the only manipulator)
    keyswitchManipulator->selectMatrixManipulator(0);  // Zero based index Value
}
void COsg::initCameraConfig(void)
{
    // Local Variable to hold window size data
    RECT rect;

    // Create the viewer for this window
    mViewer = new osgViewer::Viewer();

    // Add a Stats Handler to the viewer
    mViewer->addEventHandler(new osgViewer::StatsHandler);

	mViewer->addEventHandler(new PickHandler(mRoot.get()));  //�������
    
    // Get the current window size
    ::GetWindowRect(m_hWnd, &rect);

    // Init the GraphicsContext Traits
    osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;

    // Init the Windata Variable that holds the handle for the Window to display OSG in.
    osg::ref_ptr<osg::Referenced> windata = new osgViewer::GraphicsWindowWin32::WindowData(m_hWnd);

    // Setup the traits parameters
    traits->x = 0;
    traits->y = 0;
    traits->width = rect.right - rect.left;
    traits->height = rect.bottom - rect.top;
    traits->windowDecoration = false;
    traits->doubleBuffer = true;
    traits->sharedContext = 0;
    traits->setInheritedWindowPixelFormat = true;
    traits->inheritedWindowData = windata;

    // Create the Graphics Context
    osg::GraphicsContext* gc = osg::GraphicsContext::createGraphicsContext(traits.get());

    // Init a new Camera (Master for this View)
    osg::ref_ptr<osg::Camera> camera = new osg::Camera;

    // Assign Graphics Context to the Camera
    camera->setGraphicsContext(gc);

    // Set the viewport for the Camera
    camera->setViewport(new osg::Viewport(traits->x, traits->y, traits->width, traits->height));

    // Set projection matrix and camera attribtues
    camera->setClearMask(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
    camera->setClearColor(osg::Vec4f(0.2f, 0.2f, 0.4f, 1.0f));
    camera->setProjectionMatrixAsPerspective(
        30.0f, static_cast<double>(traits->width)/static_cast<double>(traits->height), 1.0, 1000.0);

    // Add the Camera to the Viewer
    //mViewer->addSlave(camera.get());
    mViewer->setCamera(camera.get());

    // Add the Camera Manipulator to the Viewer
    mViewer->setCameraManipulator(keyswitchManipulator.get());

    // Set the Scene Data
    mViewer->setSceneData(mRoot.get());

    // Realize the Viewer
    mViewer->realize();

    // Correct aspect ratio
    /*double fovy,aspectRatio,z1,z2;
    mViewer->getCamera()->getProjectionMatrixAsPerspective(fovy,aspectRatio,z1,z2);
    aspectRatio=double(traits->width)/double(traits->height);
    mViewer->getCamera()->setProjectionMatrixAsPerspective(fovy,aspectRatio,z1,z2);*/
}

void COsg::Render(void* ptr)
{
    COsg* osg = (COsg*)ptr;

    osgViewer::Viewer* viewer = osg->getViewer();

    // You have two options for the main viewer loop
    //      viewer->run()   or
    //      while(!viewer->done()) { viewer->frame(); }

    //viewer->run();
    while(!viewer->done())
    {
        //osg->PreFrameUpdate();
        viewer->frame();
        //osg->PostFrameUpdate();
        //Sleep(10);         // Use this command if you need to allow other processes to have cpu time
    }

    // For some reason this has to be here to avoid issue: 
    // if you have multiple OSG windows up 
    // and you exit one then all stop rendering
    //AfxMessageBox("Exit Rendering Thread");

    _endthread();
}

osg::ref_ptr<osg::Node> COsg::readPly2(std::string fileName)
{
	//��ȡ�ļ�
	FILE * pfile;
	pfile = fopen(fileName.c_str(), "r");
	int nVertex;
	int nFace;
	fscanf(pfile, "%d%d", &nVertex, &nFace);
	osg::ref_ptr<osg::Vec3Array> vertexList = new osg::Vec3Array();
	osg::ref_ptr<osg::DrawElementsUInt> faceList = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
	

	for(int i = 0; i < nVertex; i++) {
		double tmpX, tmpY, tmpZ;
		fscanf(pfile, "%lf%lf%lf", &tmpX, &tmpY, &tmpZ);
		vertexList->push_back(*(new osg::Vec3(tmpX, tmpY, tmpZ)));
	}
	for(int i = 0; i < nFace; i++) {
		int tmpX, tmpY, tmpZ, tmp;
		fscanf(pfile, "%d%d%d%d", &tmp, &tmpX, &tmpY, &tmpZ);
		faceList->push_back(tmpX);
		faceList->push_back(tmpY);
		faceList->push_back(tmpZ);
	}
	osg::ref_ptr<osg::Vec3Array> normalList = new osg::Vec3Array();  //��������Ƭ�ķ���
	normalList = calcNormalList(vertexList, faceList);
	fclose(pfile);
	//����ģ��
	osg::ref_ptr<osg::Geode> model = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	geom->setVertexArray(vertexList.get());
	geom->setNormalArray(normalList.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	geom->addPrimitiveSet(faceList.get());  //��448
	geom->setUseVertexBufferObjects(true);
	model->addDrawable(geom);
	return model;
}

osg::ref_ptr<osg::Node> COsg::readPwn(std::string pwnFileName)
{
	//��ȡ�ļ�
	FILE * pfile;
	pfile = fopen(pwnFileName.c_str(), "r");
	int nVertex;
	fscanf(pfile, "%d", &nVertex);
	osg::ref_ptr<osg::Vec3Array> vertexList = new osg::Vec3Array();	
	osg::ref_ptr<osg::Vec4Array> colorList = new osg::Vec4Array();
	osg::ref_ptr<osg::Vec3Array> normalList = new osg::Vec3Array();
	for(int i = 0; i < nVertex; i++) {
		double tmpX, tmpY, tmpZ;
		fscanf(pfile, "%lf%lf%lf", &tmpX, &tmpY, &tmpZ);
		vertexList->push_back(*(new osg::Vec3(tmpX, tmpY, tmpZ)));
		colorList->push_back(osg::Vec4(1.0, 0.0, 0.0, 1.0));
	}
	for(int i = 0; i < nVertex; i++) {
		double tmpX, tmpY, tmpZ;
		fscanf(pfile, "%lf%lf%lf", &tmpX, &tmpY, &tmpZ);
		normalList->push_back(*(new osg::Vec3(tmpX, tmpY, tmpZ)));
	}
	fclose(pfile);
	//����ģ��
	osg::ref_ptr<osg::Geode> model = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	geom->setVertexArray(vertexList.get());
	geom->setColorArray(colorList.get());
	geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
	geom->setNormalArray(normalList.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, nVertex));
	geom->setUseVertexBufferObjects(true);
	/*
	osg::Point *point=new osg::Point;
	point->setSize(0); 
	model->getOrCreateStateSet()->setAttribute(point);
	*/
	model->addDrawable(geom);
	return model;
}
osg::ref_ptr<osg::Vec3Array> COsg::calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList)
{
	//��Ҫ����ÿ������ķ���
	//�ȼ���ÿ��������Ƭ�ķ���Ȼ���ۼӣ����λ��
	osg::ref_ptr<osg::Vec3Array> vertexNormalList = new osg::Vec3Array(vertexList->size());  //����ķ���
	//�����㷨���б�����
	for(unsigned int i = 0; i < vertexList->size(); i++) {
		(*vertexNormalList)[i].set(0.0, 0.0, 0.0);
		//vertexNormalList->push_back(osg::Vec3(0.0, 0.0, 0.0));
	}
	//����ÿ����������Ƭ�ķ��򣬲��ѷ���ӵ���Ӧ����ķ���
	for(unsigned int i = 0; i < faceList->size(); i += 3) {
		osg::Vec3 vertexNormal;
		int faceIndex0, faceIndex1, faceIndex2;
		faceIndex0 = (*faceList)[i + 0];
		faceIndex1 = (*faceList)[i + 1];
		faceIndex2 = (*faceList)[i + 2];
		vertexNormal = calcFaceNormal((*vertexList)[faceIndex0], (*vertexList)[faceIndex1], (*vertexList)[faceIndex2]);
		//��Ҫ��λ��
		(*vertexNormalList)[faceIndex0] += vertexNormal;
		(*vertexNormalList)[faceIndex1] += vertexNormal;
		(*vertexNormalList)[faceIndex2] += vertexNormal;
	}
	//������ķ���λ��
	//���񲻵�λ��Ҳ���ԣ����뵥λ�����������������
	for(unsigned int i = 0; i < vertexList->size(); i++) {
		osg::Vec3 tmpVec3;
		tmpVec3 = (*vertexNormalList)[i];
		tmpVec3.normalize();
		(*vertexNormalList)[i] = tmpVec3;
	}
	return vertexNormalList;
}
osg::Vec3 COsg::calcFaceNormal(osg::Vec3 vertex0, osg::Vec3 vertex1, osg::Vec3 vertex2)
{
	osg::Vec3 tmpVector0, tmpVector1, tmpVector;
	tmpVector0 = vertex0 - vertex1;
	tmpVector1 = vertex0 - vertex2;
	tmpVector = tmpVector0 ^ tmpVector1;
	tmpVector.normalize();
	return tmpVector;
}

