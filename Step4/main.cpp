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

osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint);  //使用几何体
osg::ref_ptr<osg::Node> addDragger(osg::Node * scene);
class PickHandler : public osgGA::GUIEventHandler {
public:
	PickHandler() {}
	PickHandler(osg::Group * root) : mRoot(root) {}
	~PickHandler() {}
	bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa);
protected:
	osg::Group * mRoot;
};
int main(void)
{
	osgViewer::Viewer viewer;
	osg::ref_ptr<osg::Group> root = new osg::Group();
	osg::ref_ptr<osg::Node> tmpNode = new osg::Node();
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	viewer.addEventHandler(new PickHandler(root.get()));
	/*
	//设置点的大小
	osg::Point * point = new osg::Point();
	point->setSize(4);
	root->getOrCreateStateSet()->setAttribute(point);
	*/
	//添加点
	osg::Vec3 tmpPoint;
	tmpPoint.set(10.0, 10.0, 10.0);
	tmpNode = createPoint2(tmpPoint);
	tmpNode = addDragger(tmpNode.get());
	root->addChild(tmpNode.get());

	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	osgDB::writeNodeFile(*root, "d:\\osg2014\\step4.osg");
	//测试，root的层次关系。成功
	//std::cout << "root's NumChildren" << root->getChild(0)->asGroup()->getChild(0)->asGroup()->getNumChildren() << std::endl;

	//得到移动向量
	/*
	osg::MatrixList worldMatrices = root->getChild(0)->asGroup()->getChild(0)->getWorldMatrices();  //终于可以获得移动向量了
	for(osg::MatrixList::iterator itr = worldMatrices.begin();
		itr != worldMatrices.end();
		++itr) {
		osg::Matrix& matrix = *itr;
		osg::Vec3 trans = matrix.getTrans();
		std::cout << "move matrix: " << trans.x() << ", " << trans.y() << ", " << trans.z() << std::endl;
	}
	*/
	return 0;
}
osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint)
{
	//使用几何体
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	osg::TessellationHints * hints = new osg::TessellationHints();
	hints->setDetailRatio(0.1f);  //设置几何体的精细度，越大越精细，书上用的0.5
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
	//设置拖拽器类型
	osgManipulator::TranslateAxisDragger* d = new osgManipulator::TranslateAxisDragger();
	d->setupDefaultGeometry();
	osgManipulator::Dragger* dragger = d;
	
	//适配尺寸
	root->addChild(dragger);

	float scale = scene->getBound().radius() * 1.6;
	dragger->setMatrix(osg::Matrix::scale(scale, scale, scale) * osg::Matrix::translate(scene->getBound().center()));
	dragger->addTransformUpdating(selection);
	dragger->setHandleEvents(true);

	//在此添加激活键
	return root;
}
bool PickHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
{
    switch(ea.getEventType())
    {
        case(osgGA::GUIEventAdapter::KEYDOWN):
        {
            if (ea.getKey()=='r') {        
                //刷新并获得移动向量
				osg::MatrixList worldMatrices = mRoot->getChild(0)->asGroup()->getChild(0)->getWorldMatrices();  //终于可以获得移动向量了
				for(osg::MatrixList::iterator itr = worldMatrices.begin();
					itr != worldMatrices.end();
					++itr) {
					osg::Matrix& matrix = *itr;
					osg::Vec3 trans = matrix.getTrans();
					std::cout << "move matrix: " << trans.x() << ", " << trans.y() << ", " << trans.z() << std::endl;
				}
            }
            return false;
        }    
        default:
            return false;
    }
}