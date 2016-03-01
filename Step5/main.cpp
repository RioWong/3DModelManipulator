#include <osgViewer/Viewer>
#include <osg/Node>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/MatrixTransform>
#include <iostream>
osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint);  //创建几何体
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
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	viewer.addEventHandler(new PickHandler(root.get()));
	root->setDataVariance(osg::Object::DYNAMIC);  //场景数据的动态更新
	osg::Vec3 tmpPoint;
	tmpPoint.set(0.0, 0.0, 0.0);
	root->addChild(createPoint2(tmpPoint));
	tmpPoint.set(0.0, 4.0, 0.0);
	root->addChild(createPoint2(tmpPoint));
	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	osgDB::writeNodeFile(*root, "d:\\osg2014\\step5.osg");
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
bool PickHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
{
    switch(ea.getEventType())
    {
        case(osgGA::GUIEventAdapter::KEYDOWN):
        {
            if (ea.getKey()=='g') {        
                //给定向量，位置变动
				std::cout << "g pressed." << std::endl;
				/*
				mRoot->removeChild(mRoot->getChild(0));  //移动可以通过删除新建完成
				osg::Vec3 tmpPoint;
				tmpPoint.set(0.0, 0.0, 3.0);
				mRoot->addChild(createPoint2(tmpPoint));
				*/
				//osg::MatrixTransform * mt = dynamic_cast<osg::MatrixTransform *>(mRoot);
				
				osg::ref_ptr<osg::MatrixTransform> mt = new osg::MatrixTransform();
				osg::Matrix tmpMatrix;
				tmpMatrix.makeTranslate(-2.0, 0.0, 0.0);
				mt->setMatrix(tmpMatrix);
				mt->addChild(mRoot->getChild(1));
				mRoot->removeChild(mRoot->getChild(1));
				mRoot->addChild(mt);
            }
            return false;
        }    
        default:
            return false;
    }
}