#pragma once
#include "Osg.h"
class PickHandler : public osgGA::GUIEventHandler
{
public:
	PickHandler(void);
	~PickHandler(void);
	PickHandler(osg::Group * root);
	bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa);
	virtual void pick(osgViewer::View* view, const osgGA::GUIEventAdapter& ea);
	osg::Vec3 find(osg::Vec3 pickPoint, osg::Vec3Array * allPoint);
protected:
	osg::Vec3 resultPoint;
	osg::Group * mRoot;
private:
	double calcDistance(osg::Vec3 point1, osg::Vec3 point2);
	osg::ref_ptr<osg::Vec3Array> getNearVertexList(osg::Vec3 givenVertex, float distanceThreshold, osg::ref_ptr<osg::Vec3Array> allVertexList);
	osg::Vec3 getNormalByNear(osg::ref_ptr<osg::Vec3Array> nearVertexList);
	osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint);
	osg::ref_ptr<osg::Node> addDragger(osg::Node * scene);
};

