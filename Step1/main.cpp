#include <osgViewer/Viewer>
#include <osg/Node>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgViewer/ViewerEventHandlers>
#include <osgManipulator/TranslateAxisDragger>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/Point>
#include <iostream>

osg::ref_ptr<osg::Node> createPoint(osg::Vec3 tmpPoint);  //ֻʹ�õ�
osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint);  //ʹ�ü�����
osg::ref_ptr<osg::Node> addDragger(osg::Node * scene);
int main(void)
{
	osgViewer::Viewer viewer;
	osg::ref_ptr<osg::Group> root = new osg::Group();
	osg::ref_ptr<osg::Node> tmpNode = new osg::Node();
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	/*
	//���õ�Ĵ�С
	osg::Point * point = new osg::Point();
	point->setSize(4);
	root->getOrCreateStateSet()->setAttribute(point);
	*/
	//��ӵ�
	osg::Vec3 tmpPoint;
	tmpPoint.set(1.0, 1.0, 1.0);
	tmpNode = createPoint2(tmpPoint);
	tmpNode = addDragger(tmpNode.get());
	root->addChild(tmpNode.get());
	viewer.setSceneData(root.get());
	osgDB::writeNodeFile(*root, "d:\\osg2014\\step1.osg");
	viewer.realize();
	osgDB::writeNodeFile(*root, "d:\\osg2014\\step1rl.osg");
	viewer.run();
	osgDB::writeNodeFile(*root, "d:\\osg2014\\step1run.osg");
	return 0;
}
osg::ref_ptr<osg::Node> createPoint(osg::Vec3 tmpPoint)
{
	//ֻʹ�õ�
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	//���õ��λ��
	osg::ref_ptr<osg::Vec3Array> v = new osg::Vec3Array();
	v->push_back(tmpPoint);
	geom->setVertexArray(v.get());
	//���ö���Ĺ�����ʽ
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, 1));
	//���õ����ɫ
	osg::ref_ptr<osg::Vec4Array> c = new osg::Vec4Array();
	geom->setColorArray(c.get());
	geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
	c->push_back(osg::Vec4(1.0, 1.0, 1.0, 1.0));
	//���õ�ķ�����
	osg::ref_ptr<osg::Vec3Array> n = new osg::Vec3Array;
	geom->setNormalArray(n.get());
	geom->setNormalBinding(osg::Geometry::BIND_OVERALL);
	n->push_back(osg::Vec3(0.f, -1.f, 0.f));
	//��������������
	geode->addDrawable(geom.get());
	return geode.get();
}
osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint)
{
	//ʹ�ü�����
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	osg::TessellationHints * hints = new osg::TessellationHints();
	hints->setDetailRatio(0.1f);  //���ü�����ľ�ϸ�ȣ�Խ��Խ��ϸ�������õ�0.5
	float radius = 0.5;
	geode->addDrawable(new osg::ShapeDrawable(new osg::Box(tmpPoint, radius), hints));
	return geode.get();
}
osg::ref_ptr<osg::Node> addDragger(osg::Node * scene)
{
	scene->getOrCreateStateSet()->setMode(GL_NORMALIZE, osg::StateAttribute::ON);
	osg::MatrixTransform* selection = new osg::MatrixTransform;
	selection->addChild(scene);
	osg::Group* root = new osg::Group;
	root->addChild(selection);
	//������ק������
	//osgManipulator::Dragger* dragger = createDragger(name);
	osgManipulator::TranslateAxisDragger* d = new osgManipulator::TranslateAxisDragger();
	d->setupDefaultGeometry();
	osgManipulator::Dragger* dragger = d;
	
	//����ߴ�
	root->addChild(dragger);

	float scale = scene->getBound().radius() * 1.6;
	dragger->setMatrix(osg::Matrix::scale(scale, scale, scale) * osg::Matrix::translate(scene->getBound().center()));
	dragger->addTransformUpdating(selection);
	dragger->setHandleEvents(true);
	//�ڴ���Ӽ����
	return root;
}