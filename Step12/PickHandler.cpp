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
			if (view && mRoot->getNumChildren() == 2) pick(view,ea);  //保证拖拽器的数目
            return false;
        }
		case(osgGA::GUIEventAdapter::KEYDOWN):
        {
			if (ea.getKey()=='r' || ea.getKey() == 'R') {
                //刷新并获得移动向量
				osg::MatrixList worldMatrices = mRoot->getChild(2)->asGroup()->getChild(0)->getWorldMatrices();  //终于可以获得移动向量了
				osg::Vec3 trans;
				for(osg::MatrixList::iterator itr = worldMatrices.begin();
					itr != worldMatrices.end();
					++itr) {
					osg::Matrix& matrix = *itr;
					trans = matrix.getTrans();
					std::cout << "move matrix: " << trans.x() << ", " << trans.y() << ", " << trans.z() << std::endl;
				}
				
				//邻近顶点的位置变化
				//对mRoot的模型的minIndex顶点实行trans变换
				osg::ref_ptr<osg::Geometry> oldGeom = mRoot->getChild(0)->asGeode()->getDrawable(0)->asGeometry();
				osg::ref_ptr<osg::Vec3Array> oldVertexList = dynamic_cast<osg::Vec3Array *>(oldGeom->getVertexArray());
				osg::ref_ptr<osg::Vec3Array> oldNormalList = dynamic_cast<osg::Vec3Array *>(oldGeom->getNormalArray());
				osg::ref_ptr<osg::DrawElementsUInt> oldFaceList = dynamic_cast<osg::DrawElementsUInt *>(oldGeom->getPrimitiveSet(0));
				

				//设置新属性
				osg::ref_ptr<osg::Vec3Array> newVertexList = new osg::Vec3Array(oldVertexList->size());
				osg::ref_ptr<osg::Vec3Array> newNormalList = new osg::Vec3Array(oldNormalList->size());
				osg::ref_ptr<osg::DrawElementsUInt> newFaceList = oldFaceList;
				//对属性的处理
				//对顶点的处理
				//找出邻近的点，并改变顶点坐标
				osg::Vec3 tmpVec3(1.0, 1.0, 1.0);
				//std::cout << resultPoint.x() << ", " << resultPoint.y() << ", " << resultPoint.z() << std::endl;
				newVertexList = oldVertexList;
				float distanceThreshold = 0.1;  //顶点坐标改变的距离阀值
				float nearThreshold = 0.1;  //求法向时附近点的阀值，注意使用新顶点列表
				for(unsigned i = 0; i < oldVertexList->size(); i++) {
					float tmpDistance;
					tmpDistance = calcDistance((*oldVertexList)[i], resultPoint);
					//std::cout << "tmpDistance = " << tmpDistance << std::endl;
					if(tmpDistance <= distanceThreshold) {
						//变形算法
						(*newVertexList)[i].set((*oldVertexList)[i].x() + trans.x(), 
							                    (*oldVertexList)[i].y() + trans.y(), 
												(*oldVertexList)[i].z() + trans.z());
						//break;
					}
				}
				//newVertexList = oldVertexList;


				//对法向的处理
				//找出邻近的点，改变这些顶点的法向
				//每个法向量的改变：找出附近范围的点，任意两个点与该点组成三角形面片，求面片法向量并单位化，累加所有面片法向量并单位化
				std::cout << "process begin..." << std::endl;
				for(unsigned i = 0; i < newVertexList->size(); i++) {
					float tmpDistance = calcDistance((*oldVertexList)[i], resultPoint);
					if(tmpDistance <= distanceThreshold) {
						osg::ref_ptr<osg::Vec3Array> nearVertexList = getNearVertexList((*newVertexList)[i], nearThreshold, newVertexList);  //顶点列表的第一个点为原点
						std::cout << "done nearVertexList" << std::endl;
						(*newNormalList)[i].set(getNormalByNear(nearVertexList));
					} else {
						(*newNormalList)[i] = (*oldNormalList)[i];
					}
				}
				std::cout << "process end." << std::endl;

				//newNormalList = oldNormalList;

				//创建新模型
				osg::ref_ptr<osg::Geometry> newGeom = new osg::Geometry();
				newGeom->setVertexArray(newVertexList.get());
				newGeom->setNormalArray(newNormalList.get());
				newGeom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
				newGeom->addPrimitiveSet(newFaceList.get());  //行448
				newGeom->setUseVertexBufferObjects(true);
				//替换：删除原画布，添加几何体
				mRoot->getChild(0)->asGeode()->removeDrawable(mRoot->getChild(0)->asGeode()->getDrawable(0));
				mRoot->getChild(0)->asGeode()->addDrawable(newGeom.get());
				//mRoot->addChild(readPly2("d:\\osg2014\\model\\bdf2.ply2"));  //测试成功
				//删除拖拽器，以便进行下一步操作
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
			//打印点击的点
			std::cout << "pick vertex: " << hitr->getWorldIntersectPoint().x() << ", " 
			                               << hitr->getWorldIntersectPoint().y() << ", " 
										   << hitr->getWorldIntersectPoint().z() << std::endl;
			//在这里寻找最近的点并选中
			/*  //老方法
			//遍历所有顶点
			VertexExtractor ivea;
			mRoot->accept(ivea);
			//保存所有顶点
			osg::ref_ptr<osg::Vec3Array> allPoint = new osg::Vec3Array();
			allPoint->clear();
			allPoint->insert(allPoint->end(), ivea.extracted_verts->begin(), ivea.extracted_verts->end());
			*/
			osg::ref_ptr<osg::Vec3Array> allPoint = dynamic_cast<osg::Vec3Array *>(mRoot->getChild(1)->asGeode()->getDrawable(0)->asGeometry()->getVertexArray());
			/*
			for(unsigned i= 0; i < allPoint->size(); i++) {
				printf("%f, %f, %f\n", (*allPoint)[i].x(), (*allPoint)[i].y(), (*allPoint)[i].z());
			}
			*/ //测试成功
			osg::Vec3 pickPoint;
			//World和Local有什么区别呢
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
			//测试，输出allPoint,正常
			break;
        }
    }
}
osg::Vec3 PickHandler::find(osg::Vec3 pickPoint, osg::Vec3Array * allPoint)
{
	//计算距离，返回距离最小的点
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
	//设置拖拽器类型
	//osgManipulator::Dragger* dragger = createDragger(name);
	osgManipulator::TranslateAxisDragger* d = new osgManipulator::TranslateAxisDragger();
	d->setupDefaultGeometry();
	osgManipulator::Dragger* dragger = d;
	
	//适配尺寸
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
	//在此添加激活键
	return root;
}

osg::ref_ptr<osg::Node> PickHandler::createPoint2(osg::Vec3 tmpPoint)
{
	//使用几何体
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	osg::TessellationHints * hints = new osg::TessellationHints();
	hints->setDetailRatio(0.1f);  //设置几何体的精细度，越大越精细，书上用的0.5
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