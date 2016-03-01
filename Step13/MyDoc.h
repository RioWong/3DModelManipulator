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

#include "MeshData.h"
#include "FeatureDetector.h"
#include "Smoother.h"

class MyDoc
{
public:
	MyDoc(void);
	~MyDoc(void);
public:
	osg::ref_ptr<osg::Node> readPly2(char * fileName);
	osg::ref_ptr<osg::Vec3Array> calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList);
	osg::Vec3 calcFaceNormal(osg::Vec3 vertex0, osg::Vec3 vertex1, osg::Vec3 vertex2);

	osg::ref_ptr<osg::Node> readMesh(std::string fileName);
public:
	MeshData * mesh;
	FeatureDetector * feature_detector;
	Smoother * smoother;
};

