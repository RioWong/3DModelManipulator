#include <osgViewer/Viewer>
#include <osg/Node>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Point>
int main(void)
{
	osgViewer::Viewer viewer;
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	viewer.addEventHandler(new osgViewer::StatsHandler());
	osg::ref_ptr<osg::Group> root = new osg::Group();
	osg::ref_ptr<osg::Geode> model = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	osg::ref_ptr<osg::Vec3Array> v = new osg::Vec3Array();
	v->push_back(osg::Vec3(1, 1, 1));
	v->push_back(osg::Vec3(1, 1, 2));
	v->push_back(osg::Vec3(1, 1, 3));
	geom->setVertexArray(v.get());
	//…Ë÷√—’…´
	osg::ref_ptr<osg::Vec4Array> c = new osg::Vec4Array();
	c->push_back(osg::Vec4(1.0, 0.0, 0.0, 1.0));
	geom->setColorArray(c.get());
	geom->setColorBinding(osg::Geometry::BIND_OVERALL);
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, 3));
	osg::ref_ptr<osg::Vec3Array> n = new osg::Vec3Array();
	n->push_back(osg::Vec3(0.0, -1.0, 0.0));
	geom->setNormalArray(n.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_PRIMITIVE_SET);
	model->addDrawable(geom.get());
	root->addChild(model.get());
	osg::Point *point=new osg::Point;  
	point->setSize(40); 
	model->getOrCreateStateSet()->setAttribute(point);
	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	return 0;
}
