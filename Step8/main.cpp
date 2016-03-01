#include <osgViewer/Viewer>
#include <osg/Node>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <iostream>
osg::ref_ptr<osg::Node> readPly2(char * fileName);
osg::ref_ptr<osg::Vec3Array> calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList);
osg::Vec3 calcFaceNormal(osg::Vec3 vertex0, osg::Vec3 vertex1, osg::Vec3 vertex2);
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
	root->addChild(readPly2("d:\\osg2014\\model\\bdf2.ply2"));
	//root->addChild(osgDB::readNodeFile("d:\\osg2014\\model\\cup_2.ply"));
	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	//osgDB::writeNodeFile(*root, "d:\\osg2014\\step7ply2.osg");
	return 0;
}
osg::ref_ptr<osg::Node> readPly2(char * fileName)
{
	//读取文件
	FILE * pfile;
	pfile = fopen(fileName, "r");
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
	osg::ref_ptr<osg::Vec3Array> normalList = new osg::Vec3Array();  //三角形面片的法向
	normalList = calcNormalList(vertexList, faceList);
	fclose(pfile);
	//创建模型
	osg::ref_ptr<osg::Geode> model = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	geom->setVertexArray(vertexList.get());
	geom->setNormalArray(normalList.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	geom->addPrimitiveSet(faceList.get());  //行448
	geom->setUseVertexBufferObjects(true);
	model->addDrawable(geom);
	return model;
}
osg::ref_ptr<osg::Vec3Array> calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList)
{
	//需要的是每个顶点的法向
	//先计算每个三角面片的法向，然后累加，最后单位化
	osg::ref_ptr<osg::Vec3Array> vertexNormalList = new osg::Vec3Array(vertexList->size());  //顶点的法向
	//将顶点法向列表清零
	for(unsigned int i = 0; i < vertexList->size(); i++) {
		(*vertexNormalList)[i].set(0.0, 0.0, 0.0);
		//vertexNormalList->push_back(osg::Vec3(0.0, 0.0, 0.0));
	}
	//计算每个三角形面片的法向，并把法向加到对应顶点的法向
	for(unsigned int i = 0; i < faceList->size(); i += 3) {
		osg::Vec3 vertexNormal;
		int faceIndex0, faceIndex1, faceIndex2;
		faceIndex0 = (*faceList)[i + 0];
		faceIndex1 = (*faceList)[i + 1];
		faceIndex2 = (*faceList)[i + 2];
		vertexNormal = calcFaceNormal((*vertexList)[faceIndex0], (*vertexList)[faceIndex1], (*vertexList)[faceIndex2]);
		//需要单位化
		(*vertexNormalList)[faceIndex0] += vertexNormal;
		(*vertexNormalList)[faceIndex1] += vertexNormal;
		(*vertexNormalList)[faceIndex2] += vertexNormal;
	}
	//将顶点的法向单位化
	//好像不单位化也可以？必须单位化，否则光线有问题
	for(unsigned int i = 0; i < vertexList->size(); i++) {
		osg::Vec3 tmpVec3;
		tmpVec3 = (*vertexNormalList)[i];
		tmpVec3.normalize();
		(*vertexNormalList)[i] = tmpVec3;
	}
	return vertexNormalList;
}
osg::Vec3 calcFaceNormal(osg::Vec3 vertex0, osg::Vec3 vertex1, osg::Vec3 vertex2)
{
	osg::Vec3 tmpVector0, tmpVector1, tmpVector;
	tmpVector0 = vertex0 - vertex1;
	tmpVector1 = vertex0 - vertex2;
	tmpVector = tmpVector0 ^ tmpVector1;
	tmpVector.normalize();
	return tmpVector;
}
bool PickHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
{
    switch(ea.getEventType())
    {
        case(osgGA::GUIEventAdapter::KEYDOWN):
        {
            if (ea.getKey()=='r') {        
                //刷新并获得移动向量
				//刷新并使模型变化
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
