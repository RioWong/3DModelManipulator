#include <osgViewer/Viewer>
#include <osgDB/ReadFile>
#include <osg/Node>
#include <iostream>

int main(void)
{
	osgViewer::Viewer viewer;
	osg::ref_ptr<osg::Group> root = new osg::Group();
	root->addChild(osgDB::readNodeFile("glider.osg"));
	root->addChild(osgDB::readNodeFile("glider.osg"));
	std::cout << "Hello" << std::endl;
	std::cout << root->getNumChildren() << std::endl;
	viewer.setSceneData(root);
	viewer.realize();
	viewer.run();
	return 0;
}