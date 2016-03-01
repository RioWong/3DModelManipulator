#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osg/Node>
#include <osgViewer/ViewerEventHandlers>
#include <osg/NodeVisitor>
#include <iostream>
class VertexExtractor : public osg::NodeVisitor
{
public:
	osg::ref_ptr<osg::Vec3Array> extracted_verts;
	VertexExtractor() : osg::NodeVisitor(osg::NodeVisitor::TRAVERSE_ALL_CHILDREN)
	{
		extracted_verts = new osg::Vec3Array;
	}
	void apply( osg::Geode& geode )
	{
		for( unsigned int i=0; i < geode.getNumDrawables(); ++i )
		{
			osg::Geometry* geom = dynamic_cast<osg::Geometry*>( geode.getDrawable(i) );
			if( !geom ) {
				continue;
			}
			osg::Vec3Array* verts = dynamic_cast<osg::Vec3Array*>( geom->getVertexArray() );
			if( !verts ) {
				continue;
			}
			extracted_verts->insert( extracted_verts->end(), verts->begin(), verts->end() );
		}
	}
};
class PickHandler : public osgGA::GUIEventHandler {
public:
	PickHandler() {}
	PickHandler(osg::Group * root) : mRoot(root) {}
	~PickHandler() {}
	bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa);
	virtual void pick(osgViewer::View* view, const osgGA::GUIEventAdapter& ea);
	osg::Vec3 find(osg::Vec3 pickPoint, osg::Vec3Array * allPoint);
protected:
	osg::Group * mRoot;
};
int main(void)
{
	osgViewer::Viewer viewer;
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	osg::ref_ptr<osg::Group> root = new osg::Group();
	osg::ref_ptr<osg::Node> node = new osg::Node();
	node = osgDB::readNodeFile("glider.osg");
	root->addChild(node.get());
	viewer.addEventHandler(new PickHandler(root.get()));
	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
}
bool PickHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
{
    switch(ea.getEventType())
    {
        case(osgGA::GUIEventAdapter::PUSH):
        {
            osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
            if (view) pick(view,ea);
            return false;
        }
        default:
            return false;
    }
}
void PickHandler::pick(osgViewer::View* view, const osgGA::GUIEventAdapter& ea)
{
    osgUtil::LineSegmentIntersector::Intersections intersections;
    float x = ea.getX();
    float y = ea.getY();
    if (view->computeIntersections(x,y,intersections))
    {
        for(osgUtil::LineSegmentIntersector::Intersections::iterator hitr = intersections.begin();
            hitr != intersections.end();
            ++hitr)
        {
			std::cout << "coords vertex: " << hitr->getLocalIntersectPoint().x() << ", " 
			                               << hitr->getLocalIntersectPoint().y() << ", " 
										   << hitr->getLocalIntersectPoint().z() << std::endl;
			//在这里寻找最近的点并选中
			//遍历所有顶点
			VertexExtractor ivea;
			mRoot->accept(ivea);
			//保存所有顶点
			osg::ref_ptr<osg::Vec3Array> allPoint = new osg::Vec3Array();
			allPoint->clear();
			allPoint->insert(allPoint->end(), ivea.extracted_verts->begin(), ivea.extracted_verts->end());
			osg::Vec3 pickPoint, resultPoint;
			pickPoint.set(hitr->getLocalIntersectPoint().x(), hitr->getLocalIntersectPoint().y(), hitr->getLocalIntersectPoint().z());
			resultPoint = this->find(pickPoint, allPoint);
			std::cout << "result point: " << resultPoint.x() << ", "
				                          << resultPoint.y() << ", "
										  << resultPoint.z() << std::endl;
			break;
        }
    }
}
osg::Vec3 PickHandler::find(osg::Vec3 pickPoint, osg::Vec3Array * allPoint)
{
	//计算距离，返回距离最小的点
	std::vector<osg::Vec3>::iterator iter = allPoint->begin();
	float minDistance = sqrt((iter->x() - pickPoint.x()) * (iter->x() - pickPoint.x()) + 
			                 (iter->x() - pickPoint.y()) * (iter->x() - pickPoint.y()) + 
							 (iter->x() - pickPoint.z()) * (iter->x() - pickPoint.z()));
	int minIndex = 0;
	osg::Vec3 minPoint;
	minPoint.set(iter->x(), iter->y(), iter->z());
	for(unsigned int i = 0; i < allPoint->size(); i++) {
		float tmpDistance = sqrt((iter->x() - pickPoint.x()) * (iter->x() - pickPoint.x()) + 
			                     (iter->x() - pickPoint.y()) * (iter->x() - pickPoint.y()) + 
								 (iter->x() - pickPoint.z()) * (iter->x() - pickPoint.z()));
		if(tmpDistance < minDistance) {
			minDistance = tmpDistance;
			minIndex = i;
			minPoint.set(iter->x(), iter->y(), iter->z());
		}
		iter++;
	}
	return minPoint;
}
