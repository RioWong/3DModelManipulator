#include <osgViewer/Viewer>
#include <osg/Node>
#include <osgViewer/ViewerEventHandlers>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <iostream>
#include <osgManipulator/TranslateAxisDragger>
#include <osg/Geometry>
#include <osg/ShapeDrawable>
#include <osg/NodeVisitor>
#include <osg/Point>
unsigned modelIndex = 0;
unsigned rrIndex = 1;
unsigned draggerIndex = 2;
osg::ref_ptr<osg::Node> readPly2(char * fileName);
osg::ref_ptr<osg::Node> readRr(char * rrFileName);
osg::ref_ptr<osg::Node> showRr(char * rrFileName);
osg::ref_ptr<osg::Node> readPwn(char * pwnFileName);
osg::ref_ptr<osg::Vec3Array> calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList);
osg::Vec3 calcFaceNormal(osg::Vec3 vertex0, osg::Vec3 vertex1, osg::Vec3 vertex2);
osg::ref_ptr<osg::Node> addDragger(osg::Node * scene);
osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint);
double calcDistance(osg::Vec3 point1, osg::Vec3 point2);
osg::ref_ptr<osg::Vec3Array> getNearVertexList(osg::Vec3 givenVertex, float distanceThreshold, osg::ref_ptr<osg::Vec3Array> allVertexList);
osg::Vec3 getNormalByNear(osg::ref_ptr<osg::Vec3Array> nearVertexList);
class PickHandler : public osgGA::GUIEventHandler {
public:
	PickHandler() {}
	PickHandler(osg::Group * root) : mRoot(root) {}
	~PickHandler() {}
	bool handle(const osgGA::GUIEventAdapter& ea,osgGA::GUIActionAdapter& aa);
	virtual void pick(osgViewer::View* view, const osgGA::GUIEventAdapter& ea);
	osg::Vec3 find(osg::Vec3 pickPoint, osg::Vec3Array * allPoint);
protected:
	osg::Vec3 resultPoint;
	osg::Group * mRoot;
};
int main(int argc, char ** argv)
{
	osgViewer::Viewer viewer;
	osg::ref_ptr<osg::Group> root = new osg::Group();
	//root->setDataVariance(osg::Object::DYNAMIC);
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	viewer.addEventHandler(new osgViewer::StatsHandler());
	if(argc == 0) {
		std::cout << "EROOR." << std::endl;
		return 1;
	} 
	//char * fileName = argv[1];
	
	char * ply2FileName = "bdf2.ply2";  //模型文件
	char * rrFileName = "bdf2.rr";  //谷脊点文件
	osg::ref_ptr<osg::Node> ply2Model = readPly2(ply2FileName);
	osg::ref_ptr<osg::Node> rrModel = readPwn(rrFileName);
	root->addChild(ply2Model);
	root->addChild(rrModel);
	modelIndex = root->getChildIndex(ply2Model);
	rrIndex = root->getChildIndex(rrModel);
	viewer.addEventHandler(new PickHandler(root.get()));
	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	std::cout << "saved begin..." << std::endl;
	//osgDB::writeNodeFile(*root, "bdf2.osg");
	//osgDB::writeNodeFile(*root, "bdf2.obj");
	std::cout << "saved successful." << std::endl;
	return 0;
}
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

// The DraggerContainer node is used to fix the dragger's size on the screen
class DraggerContainer : public osg::Group
{
public:
    DraggerContainer() : _draggerSize(240.0f), _active(true) {}
    
    DraggerContainer( const DraggerContainer& copy, const osg::CopyOp& copyop=osg::CopyOp::SHALLOW_COPY )
    :   osg::Group(copy, copyop),
        _dragger(copy._dragger), _draggerSize(copy._draggerSize), _active(copy._active)
    {}
    
    META_Node( osgManipulator, DraggerContainer );
    
    void setDragger( osgManipulator::Dragger* dragger )
    {
        _dragger = dragger;
        if ( !containsNode(dragger) ) addChild( dragger );
    }
    
    osgManipulator::Dragger* getDragger() { return _dragger.get(); }
    const osgManipulator::Dragger* getDragger() const { return _dragger.get(); }
    
    void setDraggerSize( float size ) { _draggerSize = size; }
    float getDraggerSize() const { return _draggerSize; }
    
    void setActive( bool b ) { _active = b; }
    bool getActive() const { return _active; }
    
    void traverse( osg::NodeVisitor& nv )
    {
        if ( _dragger.valid() )
        {
            if ( _active && nv.getVisitorType()==osg::NodeVisitor::CULL_VISITOR )
            {
                osgUtil::CullVisitor* cv = static_cast<osgUtil::CullVisitor*>(&nv);
                
                float pixelSize = cv->pixelSize(_dragger->getBound().center(), 0.48f);
                if ( pixelSize!=_draggerSize )
                {
                    float pixelScale = pixelSize>0.0f ? _draggerSize/pixelSize : 1.0f;
                    osg::Vec3d scaleFactor(pixelScale, pixelScale, pixelScale);
                    
                    osg::Vec3 trans = _dragger->getMatrix().getTrans();
                    _dragger->setMatrix( osg::Matrix::scale(scaleFactor) * osg::Matrix::translate(trans) );
                }
            }
        }
        osg::Group::traverse(nv);
    }
    
protected:
    osg::ref_ptr<osgManipulator::Dragger> _dragger;
    float _draggerSize;
    bool _active;
};

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
osg::ref_ptr<osg::Node> readRr(char * rrFileName)
{
	//只存储信息，不显示
	FILE * pfile;
	osg::ref_ptr<osg::Vec3Array> rrList = new osg::Vec3Array();
	int nRidge, nRavine;
	pfile = fopen(rrFileName, "r");
	fscanf(pfile, "%d", &nRidge);
	for(int i = 0; i < nRidge; i++) {
		osg::Vec3 tmpVec3;
		fscanf(pfile, "%f%f%f", &tmpVec3.x(), &tmpVec3.y(), &tmpVec3.z());
		rrList->push_back(tmpVec3);
		//printf("%f, %f, %f\n", tmpVec3.x(), tmpVec3.y(), tmpVec3.z());  //测试成功
	}
	fscanf(pfile, "%d", &nRavine);
	for(int i = 0; i < nRavine; i++) {
		osg::Vec3 tmpVec3;
		fscanf(pfile, "%f%f%f", &tmpVec3.x(), &tmpVec3.y(), &tmpVec3.z());
		rrList->push_back(tmpVec3);
	}
	fclose(pfile);

	osg::ref_ptr<osg::Geode> rr = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	geom->setVertexArray(rrList.get());
	rr->addDrawable(geom);
	return rr;
}
osg::ref_ptr<osg::Node> showRr(char * rrFileName)
{
	//显示谷脊点
	FILE * pfile;
	osg::ref_ptr<osg::Vec3Array> rrList = new osg::Vec3Array();
	int nRidge, nRavine;
	pfile = fopen(rrFileName, "r");
	fscanf(pfile, "%d", &nRidge);
	for(int i = 0; i < nRidge; i++) {
		osg::Vec3 tmpVec3;
		fscanf(pfile, "%f%f%f", &tmpVec3.x(), &tmpVec3.y(), &tmpVec3.z());
		rrList->push_back(tmpVec3);
		//printf("%f, %f, %f\n", tmpVec3.x(), tmpVec3.y(), tmpVec3.z());  //测试成功
	}
	fscanf(pfile, "%d", &nRavine);
	for(int i = 0; i < nRavine; i++) {
		osg::Vec3 tmpVec3;
		fscanf(pfile, "%f%f%f", &tmpVec3.x(), &tmpVec3.y(), &tmpVec3.z());
		rrList->push_back(tmpVec3);
	}
	fclose(pfile);

	osg::ref_ptr<osg::Geode> rr = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	geom->setVertexArray(rrList.get());
	//设置颜色
	osg::ref_ptr<osg::Vec4Array> c = new osg::Vec4Array();
	c->push_back(osg::Vec4(1.0, 0.0, 0.0, 1.0));
	geom->setColorArray(c.get());
	geom->setColorBinding(osg::Geometry::BIND_OVERALL);
	/*
	//定义法线
	osg::ref_ptr<osg::Vec3Array> n = new osg::Vec3Array();
	n->push_back(osg::Vec3(0.0, 1.0, 0.0));
	geom->setNormalArray(n.get());
	geom->setNormalBinding(osg::Geometry::BIND_OVERALL);
	*/
	rr->addDrawable(geom);
	osg::Point *point=new osg::Point;  
	point->setSize(40); 
	rr->getOrCreateStateSet()->setAttribute(point);
	return rr;
}
osg::ref_ptr<osg::Node> readPwn(char * pwnFileName)
{
	//读取文件
	FILE * pfile;
	pfile = fopen(pwnFileName, "r");
	int nVertex;
	fscanf(pfile, "%d", &nVertex);
	osg::ref_ptr<osg::Vec3Array> vertexList = new osg::Vec3Array();	
	osg::ref_ptr<osg::Vec4Array> colorList = new osg::Vec4Array();
	osg::ref_ptr<osg::Vec3Array> normalList = new osg::Vec3Array();
	for(int i = 0; i < nVertex; i++) {
		double tmpX, tmpY, tmpZ;
		fscanf(pfile, "%lf%lf%lf", &tmpX, &tmpY, &tmpZ);
		vertexList->push_back(*(new osg::Vec3(tmpX, tmpY, tmpZ)));
		colorList->push_back(osg::Vec4(1.0, 0.0, 0.0, 1.0));
	}
	for(int i = 0; i < nVertex; i++) {
		double tmpX, tmpY, tmpZ;
		fscanf(pfile, "%lf%lf%lf", &tmpX, &tmpY, &tmpZ);
		normalList->push_back(*(new osg::Vec3(tmpX, tmpY, tmpZ)));
	}
	fclose(pfile);
	//创建模型
	osg::ref_ptr<osg::Geode> model = new osg::Geode();
	osg::ref_ptr<osg::Geometry> geom = new osg::Geometry();
	geom->setVertexArray(vertexList.get());
	geom->setColorArray(colorList.get());
	geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);
	geom->setNormalArray(normalList.get());
	geom->setNormalBinding(osg::Geometry::BIND_PER_VERTEX);
	
	geom->addPrimitiveSet(new osg::DrawArrays(osg::PrimitiveSet::POINTS, 0, nVertex));
	geom->setUseVertexBufferObjects(true);
	/*
	osg::Point *point=new osg::Point;
	point->setSize(0); 
	model->getOrCreateStateSet()->setAttribute(point);
	*/
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
        case(osgGA::GUIEventAdapter::PUSH):
        {
            osgViewer::View* view = dynamic_cast<osgViewer::View*>(&aa);
			if (view && mRoot->getNumChildren() == 2) pick(view,ea);  //保证拖拽器的数目
            return false;
        }
		case(osgGA::GUIEventAdapter::KEYDOWN):
        {
            if (ea.getKey()=='r') {        
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
				float distanceThreshold = 1;  //顶点坐标改变的距离阀值
				float nearThreshold = 0.5;  //求法向时附近点的阀值，注意使用新顶点列表
				for(unsigned i = 0; i < oldVertexList->size(); i++) {
					float tmpDistance;
					/*
					tmpDistance = ((*oldVertexList)[i].x() - resultPoint.x()) * ((*oldVertexList)[i].x() - resultPoint.x()) +
						          ((*oldVertexList)[i].y() - resultPoint.y()) * ((*oldVertexList)[i].y() - resultPoint.y()) +
								  ((*oldVertexList)[i].z() - resultPoint.z()) * ((*oldVertexList)[i].z() - resultPoint.z());
					*/
					tmpDistance = calcDistance((*oldVertexList)[i], resultPoint);
					//std::cout << "tmpDistance = " << tmpDistance << std::endl;
					if(tmpDistance <= distanceThreshold) {
						//变形算法
						osg::Vec3 tmpVertex;
						osg::Vec3 tmpMove;
						tmpVertex.x() = (*oldVertexList)[i].x() + trans.x() * tmpDistance;
						tmpVertex.y() = (*oldVertexList)[i].y() + trans.y();
						tmpVertex.z() = (*oldVertexList)[i].z() + trans.z();
						/*
						(*newVertexList)[i].set((*oldVertexList)[i].x() + trans.x(), 
							                    (*oldVertexList)[i].y() + trans.y(), 
												(*oldVertexList)[i].z() + trans.z());
						*/
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
			/*
			std::cout << "coords vertex: " << hitr->getLocalIntersectPoint().x() << ", " 
			                               << hitr->getLocalIntersectPoint().y() << ", " 
										   << hitr->getLocalIntersectPoint().z() << std::endl;
			*/

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
osg::ref_ptr<osg::Node> addDragger(osg::Node * scene)
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
osg::ref_ptr<osg::Node> createPoint2(osg::Vec3 tmpPoint)
{
	//使用几何体
	osg::ref_ptr<osg::Geode> geode = new osg::Geode();
	osg::TessellationHints * hints = new osg::TessellationHints();
	hints->setDetailRatio(0.1f);  //设置几何体的精细度，越大越精细，书上用的0.5
	float radius = 0.3;
	geode->addDrawable(new osg::ShapeDrawable(new osg::Box(tmpPoint, radius), hints));
	return geode.get();
}
double calcDistance(osg::Vec3 point1, osg::Vec3 point2)
{
	double d;
	d = (point1.x() - point2.x()) * (point1.x() - point2.x()) +
		(point1.y() - point2.y()) * (point1.y() - point2.y()) +
		(point1.z() - point2.z()) * (point1.z() - point2.z());
	return sqrt(d);
}
osg::ref_ptr<osg::Vec3Array> getNearVertexList(osg::Vec3 givenVertex, float nearThreshold, osg::ref_ptr<osg::Vec3Array> allVertexList)
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
osg::Vec3 getNormalByNear(osg::ref_ptr<osg::Vec3Array> nearVertexList)
{
	osg::Vec3 normal;
	normal.set(0.0, 0.0, 0.0);
	osg::Vec3 givenVertex = (*nearVertexList)[0];
	for(unsigned i = 1; i < nearVertexList->size(); i++) {
		for(unsigned j = i + 1; j < nearVertexList->size(); j++) {
			osg::Vec3 faceNormal = calcFaceNormal(givenVertex, (*nearVertexList)[i], (*nearVertexList)[j]);
			faceNormal.normalize();
			normal += faceNormal;
		}
	}
	normal.normalize();
	return normal;
}
