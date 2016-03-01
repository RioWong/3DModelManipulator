#include "StdAfx.h"
#include "PickHandler.h"
#include "MainFrm.h"

PickHandler::PickHandler(void)
{
}


PickHandler::~PickHandler(void)
{
}

PickHandler::PickHandler(osg::Group * root) : mRoot(root) {}

bool PickHandler::handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa)
{
    switch(ea.getEventType())
    {
        case(osgGA::GUIEventAdapter::PUSH):
        {
			osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
			if (view && mRoot->getNumChildren() == 2) pick(view,ea);  //��֤��ק������Ŀ
            return false;
        }
		case(osgGA::GUIEventAdapter::KEYDOWN):
        {
			if (ea.getKey()=='r' || ea.getKey() == 'R') {
                //ˢ�²�����ƶ�����
				osg::MatrixList worldMatrices = mRoot->getChild(2)->asGroup()->getChild(0)->getWorldMatrices();  //���ڿ��Ի���ƶ�������
				osg::Vec3 trans;
				for(osg::MatrixList::iterator itr = worldMatrices.begin();
					itr != worldMatrices.end();
					++itr) {
					osg::Matrix& matrix = *itr;
					trans = matrix.getTrans();
					std::cout << "move matrix: " << trans.x() << ", " << trans.y() << ", " << trans.z() << std::endl;
				}
				
				//�ڽ������λ�ñ仯
				//��mRoot��ģ�͵�minIndex����ʵ��trans�任
				osg::ref_ptr<osg::Geometry> oldGeom = mRoot->getChild(0)->asGeode()->getDrawable(0)->asGeometry();
				osg::ref_ptr<osg::Vec3Array> oldVertexList = dynamic_cast<osg::Vec3Array *>(oldGeom->getVertexArray());
				osg::ref_ptr<osg::Vec3Array> oldNormalList = dynamic_cast<osg::Vec3Array *>(oldGeom->getNormalArray());
				osg::ref_ptr<osg::DrawElementsUInt> oldFaceList = dynamic_cast<osg::DrawElementsUInt *>(oldGeom->getPrimitiveSet(0));
				

				//����������
				osg::ref_ptr<osg::Vec3Array> newVertexList = new osg::Vec3Array(oldVertexList->size());
				osg::ref_ptr<osg::Vec3Array> newNormalList = new osg::Vec3Array(oldNormalList->size());
				osg::ref_ptr<osg::DrawElementsUInt> newFaceList = oldFaceList;
				//�����ԵĴ���
				//�Զ���Ĵ���
				//�ҳ��ڽ��ĵ㣬���ı䶥������
				osg::Vec3 tmpVec3(1.0, 1.0, 1.0);
				//std::cout << resultPoint.x() << ", " << resultPoint.y() << ", " << resultPoint.z() << std::endl;
				newVertexList = oldVertexList;
				float distanceThreshold = 0.1;  //��������ı�ľ��뷧ֵ
				float nearThreshold = 0.1;  //����ʱ������ķ�ֵ��ע��ʹ���¶����б�
				for(unsigned i = 0; i < oldVertexList->size(); i++) {
					float tmpDistance;
					tmpDistance = calcDistance((*oldVertexList)[i], resultPoint);
					//std::cout << "tmpDistance = " << tmpDistance << std::endl;
					if(tmpDistance <= distanceThreshold) {
						//�����㷨
						(*newVertexList)[i].set((*oldVertexList)[i].x() + trans.x(), 
							                    (*oldVertexList)[i].y() + trans.y(), 
												(*oldVertexList)[i].z() + trans.z());
						//break;
					}
				}
				//newVertexList = oldVertexList;


				//�Է���Ĵ���
				//�ҳ��ڽ��ĵ㣬�ı���Щ����ķ���
				//ÿ���������ĸı䣺�ҳ�������Χ�ĵ㣬������������õ������������Ƭ������Ƭ����������λ�����ۼ�������Ƭ����������λ��
				std::cout << "process begin..." << std::endl;
				for(unsigned i = 0; i < newVertexList->size(); i++) {
					float tmpDistance = calcDistance((*oldVertexList)[i], resultPoint);
					if(tmpDistance <= distanceThreshold) {
						osg::ref_ptr<osg::Vec3Array> nearVertexList = getNearVertexList((*newVertexList)[i], nearThreshold, newVertexList);  //�����б�ĵ�һ����Ϊԭ��
						std::cout << "done nearVertexList" << std::endl;
						(*newNormalList)[i].set(getNormalByNear(nearVertexList));
					} else {
						(*newNormalList)[i] = (*oldNormalList)[i];
					}
				}
				std::cout << "process end." << std::endl;

				//newNormalList = oldNormalList;

				//������ģ��
				osg::ref_ptr<osg::Geometry> newGeom = new osg::Geometry();
				newGeom->setVertexArray(newVertexList.get());
				newGeom->setNormalArray(newNormalList.get());
				newGeom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
				newGeom->addPrimitiveSet(newFaceList.get());  //��448
				newGeom->setUseVertexBufferObjects(true);
				//�滻��ɾ��ԭ��������Ӽ�����
				mRoot->getChild(0)->asGeode()->removeDrawable(mRoot->getChild(0)->asGeode()->getDrawable(0));
				mRoot->getChild(0)->asGeode()->addDrawable(newGeom.get());
				//mRoot->addChild(readPly2("d:\\osg2014\\model\\bdf2.ply2"));  //���Գɹ�
				//ɾ����ק�����Ա������һ������
				mRoot->removeChild(mRoot->getChild(2));
				std::cout << "manipulate successful." << std::endl;
				std::cout << "----------" << std::endl;
            }
            return false;
        }
        default:
            return false;
    }
	return false;
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
			//��ӡ����ĵ�
			std::cout << "pick vertex: " << hitr->getWorldIntersectPoint().x() << ", " 
			                               << hitr->getWorldIntersectPoint().y() << ", " 
										   << hitr->getWorldIntersectPoint().z() << std::endl;
			//������Ѱ������ĵ㲢ѡ��
			/*  //�Ϸ���
			//�������ж���
			VertexExtractor ivea;
			mRoot->accept(ivea);
			//�������ж���
			osg::ref_ptr<osg::Vec3Array> allPoint = new osg::Vec3Array();
			allPoint->clear();
			allPoint->insert(allPoint->end(), ivea.extracted_verts->begin(), ivea.extracted_verts->end());
			*/
			osg::ref_ptr<osg::Vec3Array> allPoint = dynamic_cast<osg::Vec3Array *>(mRoot->getChild(1)->asGeode()->getDrawable(0)->asGeometry()->getVertexArray());
			/*
			for(unsigned i= 0; i < allPoint->size(); i++) {
				printf("%f, %f, %f\n", (*allPoint)[i].x(), (*allPoint)[i].y(), (*allPoint)[i].z());
			}
			*/ //���Գɹ�
			osg::Vec3 pickPoint;
			//World��Local��ʲô������
			pickPoint.set(hitr->getWorldIntersectPoint().x(), hitr->getWorldIntersectPoint().y(), hitr->getWorldIntersectPoint().z());
			resultPoint = this->find(pickPoint, allPoint);
			std::cout << "result point: " << resultPoint.x() << ", "
				                          << resultPoint.y() << ", "
										  << resultPoint.z() << std::endl;
			std::cout << "get action vertex successful." << std::endl;
			std::cout << "----------" << std::endl;
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
osg::ref_ptr<osg::Node> PickHandler::addDragger(osg::Node * scene)
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
	if(true) {
		DraggerContainer* draggerContainer = new DraggerContainer;
        draggerContainer->setDragger( dragger );
        root->addChild(draggerContainer);
	} else {
		root->addChild(dragger);
	}
	
	float scale = scene->getBound().radius() * 1.6;
	dragger->setMatrix(osg::Matrix::scale(scale, scale, scale) * osg::Matrix::translate(scene->getBound().center()));
	dragger->addTransformUpdating(selection);
	dragger->setHandleEvents(true);
	//�ڴ���Ӽ����
	return root;
}

osg::ref_ptr<osg::Node> PickHandler::createPoint2(osg::Vec3 tmpPoint)
{
	//ʹ�ü�����
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	osg::TessellationHints * hints = new osg::TessellationHints();
	hints->setDetailRatio(0.1f);  //���ü�����ľ�ϸ�ȣ�Խ��Խ��ϸ�������õ�0.5
	float radius = 0.3;
	geode->addDrawable(new osg::ShapeDrawable(new osg::Box(tmpPoint, radius), hints));
	return geode.get();
}
double PickHandler::calcDistance(osg::Vec3 point1, osg::Vec3 point2)
{
	double d;
	d = (point1.x() - point2.x()) * (point1.x() - point2.x()) +
		(point1.y() - point2.y()) * (point1.y() - point2.y()) +
		(point1.z() - point2.z()) * (point1.z() - point2.z());
	return sqrt(d);
}
osg::ref_ptr<osg::Vec3Array> PickHandler::getNearVertexList(osg::Vec3 givenVertex, float nearThreshold, osg::ref_ptr<osg::Vec3Array> allVertexList)
{
	osg::ref_ptr<osg::Vec3Array> nearVertexList = new osg::Vec3Array();
	nearVertexList->push_back(givenVertex);
	for(unsigned i = 0; i < allVertexList->size(); i++) {
		if(calcDistance(givenVertex, (*allVertexList)[i]) <= nearThreshold) {
			nearVertexList->push_back((*allVertexList)[i]);
		}
	}
	return nearVertexList;
}
osg::Vec3 PickHandler::getNormalByNear(osg::ref_ptr<osg::Vec3Array> nearVertexList)
{
	osg::Vec3 normal;
	normal.set(0.0, 0.0, 0.0);
	osg::Vec3 givenVertex = (*nearVertexList)[0];
	for(unsigned i = 1; i < nearVertexList->size(); i++) {
		for(unsigned j = i + 1; j < nearVertexList->size(); j++) {
			osg::Vec3 faceNormal = COsg::calcFaceNormal(givenVertex, (*nearVertexList)[i], (*nearVertexList)[j]);
			faceNormal.normalize();
			normal += faceNormal;
		}
	}
	normal.normalize();
	return normal;
}