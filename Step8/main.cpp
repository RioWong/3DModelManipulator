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
	//��ȡ�ļ�
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
osg::ref_ptr<osg::Vec3Array> calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList)
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
                //ˢ�²�����ƶ�����
				//ˢ�²�ʹģ�ͱ仯
				osg::MatrixList worldMatrices = mRoot->getChild(0)->asGroup()->getChild(0)->getWorldMatrices();  //���ڿ��Ի���ƶ�������
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
