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

osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint);  //ʹ�ü�����
osg::ref_ptr<osg::Node> addDragger(osg::Node * scene);
osg::ref_ptr<osg::Group> readRr(char * fileName);
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
	char * fileName = "d:\\osg2014\\model\\cup_2.rr";
	osg::ref_ptr<osg::Group> root = readRr(fileName);

	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	viewer.addEventHandler(new PickHandler(root.get()));

	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	//osgDB::writeNodeFile(*root, "d:\\osg2014\\step6.osg");

	return 0;
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
        case(osgGA::GUIEventAdapter::KEYDOWN):
        {
            if (ea.getKey()=='r') {        
                //ˢ�²�����ƶ�����
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
			//������Ѱ������ĵ㲢ѡ��
			//�������ж���
			VertexExtractor ivea;
			mRoot->accept(ivea);
			//�������ж���
			osg::ref_ptr<osg::Vec3Array> allPoint = new osg::Vec3Array();
			allPoint->clear();
			allPoint->insert(allPoint->end(), ivea.extracted_verts->begin(), ivea.extracted_verts->end());
			osg::Vec3 pickPoint, resultPoint;
			pickPoint.set(hitr->getLocalIntersectPoint().x(), hitr->getLocalIntersectPoint().y(), hitr->getLocalIntersectPoint().z());
			resultPoint = this->find(pickPoint, allPoint);
			std::cout << "result point: " << resultPoint.x() << ", "
				                          << resultPoint.y() << ", "
										  << resultPoint.z() << std::endl;
			osg::ref_ptr<osg::Node> tmpNode = new osg::Node();
			tmpNode = createPoint2(resultPoint);
			//tmpNode = createPoint2(pickPoint);
			mRoot->addChild(addDragger(tmpNode));
			//���ԣ����allPoint,����
			break;
        }
    }
}
osg::Vec3 PickHandler::find(osg::Vec3 pickPoint, osg::Vec3Array * allPoint)
{
	//������룬���ؾ�����С�ĵ�
	std::vector<osg::Vec3>::iterator iter = allPoint->begin();
	double minDistance =    ((iter->x() - pickPoint.x()) * (iter->x() - pickPoint.x()) + 
			                 (iter->y() - pickPoint.y()) * (iter->y() - pickPoint.y()) + 
							 (iter->z() - pickPoint.z()) * (iter->z() - pickPoint.z()));
	int minIndex = 0;
	osg::Vec3 minPoint;
	minPoint.set(iter->x(), iter->y(), iter->z());
	for(unsigned int i = 0; i < allPoint->size(); i++) {
		double tmpDistance =    ((iter->x() - pickPoint.x()) * (iter->x() - pickPoint.x()) + 
			                     (iter->y() - pickPoint.y()) * (iter->y() - pickPoint.y()) + 
								 (iter->z() - pickPoint.z()) * (iter->z() - pickPoint.z()));
		if(tmpDistance < minDistance) {
			minDistance = tmpDistance;
			minIndex = i;
			minPoint.set(iter->x(), iter->y(), iter->z());
		}
		iter++;
	}
	return minPoint;
}
osg::ref_ptr<osg::Group> readRr(char * fileName)
{
	FILE * pfile;
	pfile = fopen(fileName, "r");
	int nRidge;
	fscanf(pfile, "%d", &nRidge);
	osg::ref_ptr<osg::Vec3Array> ridgeList = new osg::Vec3Array(nRidge);
	for(int i = 0; i < nRidge; i++) {
		float tmpX, tmpY, tmpZ;
		fscanf(pfile, "%f%f%f", &tmpX, &tmpY, &tmpZ);
		(*ridgeList)[i].set(tmpX, tmpY, tmpZ);
	}
	fclose(pfile);

	osg::ref_ptr<osg::Group> rrModel = new osg::Group();
	for(unsigned int i = 0; i < ridgeList->size(); i++) {
		rrModel->addChild(createPoint2((*ridgeList)[i]));
	}

	return rrModel;
}