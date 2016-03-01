#include "MyDoc.h"


MyDoc::MyDoc(void)
{
	/*
	osgViewer::Viewer viewer;
	viewer.addEventHandler(new osgViewer::WindowSizeHandler());
	osg::ref_ptr<osg::Group> root = new osg::Group();
	osg::ref_ptr<osg::Node> ply2Model = readPly2("d:\\osg2014\\model\\bdf2.ply2");
	root->addChild(ply2Model.get());
	viewer.setSceneData(root.get());
	viewer.realize();
	viewer.run();
	*/
	readMesh("d:\\osg2014\\model\\bdf2.ply2");
}

MyDoc::~MyDoc(void)
{
}

osg::ref_ptr<osg::Node> MyDoc::readPly2(char * fileName)
{
	//读取文件
	//还要存入MeshData
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

osg::ref_ptr<osg::Vec3Array> MyDoc::calcNormalList(osg::ref_ptr<osg::Vec3Array> vertexList, osg::ref_ptr<osg::DrawElementsUInt> faceList)
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

osg::Vec3 MyDoc::calcFaceNormal(osg::Vec3 vertex0, osg::Vec3 vertex1, osg::Vec3 vertex2)
{
	osg::Vec3 tmpVector0, tmpVector1, tmpVector;
	tmpVector0 = vertex0 - vertex1;
	tmpVector1 = vertex0 - vertex2;
	tmpVector = tmpVector0 ^ tmpVector1;
	tmpVector.normalize();
	return tmpVector;
}

osg::ref_ptr<osg::Node> MyDoc::readMesh(std::string fileName)
{
	//文件打开开始
	FILE * file;
	file = fopen(fileName.c_str(), "r");

	int vertex_N, face_N;
	fscanf(file, "%d", &vertex_N);
	mesh->setVertexCount(vertex_N);
	fscanf(file, "%d", &face_N);
	mesh->setFaceCount(face_N);

	float p[3];
	float dummy;
	int i;  for(i = 0; i<vertex_N; i++){
		fscanf(file, "%f", &p[0]);
		fscanf(file, "%f", &p[1]);
		fscanf(file, "%f", &p[2]);
		mesh->setVertex(i, p);
	}

	int f[3];
	int i;  for(i = 0; i<face_N; i++){
		fscanf(file, "%f", &dummy);
		fscanf(file, "%d", &f[0]);
		fscanf(file, "%d", &f[1]);
		fscanf(file, "%d", &f[2]);
		mesh->setFace(i, f);
	}
	//smoother->setMeshData(mesh);
	feature_detector->setMeshData(mesh);
	fclose(file);
	smoother->setMeshData(mesh);
	feature_detector->setMeshData(mesh);
	mesh->computeFaceNormal();  //计算面的法向
	mesh->computeNormal();  //计算点的法向
	//文件打开结束
	//提取谷脊线开始
	MeshData* current_mesh = NULL;
	mesh->computePrincipal();  //计算顶点主曲率最大值、最小值、方向
	feature_detector->ridge();
	//feature_detector->automaticThresholding(ridge_dialog->m_hi,ridge_dialog->m_low);
	feature_detector->automaticThresholding(0.6, 0.3);
	feature_detector->crossThresholding();  //Abosolute value comparison
	feature_detector->addBridgeRR(0, 0);  //不知何用
	//基本操作结束
	//张量平滑开始
	/*
	int iteration = 10;
	int in_iter = 1;
	int out_iter = iteration/in_iter;
	int mod_iter = iteration - out_iter*in_iter;
	int i;  for(i = 0; i<out_iter; i++){
		smoother->smoothKmaxKmin2(in_iter, 0.1);
	}
	smoother->smoothKmaxKmin2(mod_iter, 0.1);

	iteration = 5;
	in_iter = 1;
	out_iter = iteration/in_iter;
	mod_iter = iteration - out_iter*in_iter;
	int i;  for(i = 0; i<out_iter; i++){
		smoother->smoothTmaxTmin3(in_iter, 5);
	}
	smoother->smoothTmaxTmin3(in_iter, 5);
	*/
	//张量平滑结束
	//保存谷脊点数目及坐标开始
	int f_N = mesh->countValidFace();
	int v_N = mesh->countValidVertex();
	int* vertex_table = new int[mesh->vertex_N];
	int counter = 0;

	float (*vertex)[3] = mesh->vertex;
	int (*face)[3] = mesh->face;
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;

	int i;  for(i = 0; i<vertex_N; i++){
		if(mesh->degree_v[i] != 0){
			vertex_table[i] = counter;
			counter++;
		}
		else
			vertex_table[i] = -1;
	}

	int ridge_N = 0;
	int ravine_N = 0;
	int i;  for(i = 0; i<vertex_N; i++){
		if(mesh->isRidge[i])
			ridge_N++;
		if(mesh->isRavine[i])
			ravine_N++;
	}
	fprintf(file, "%d\n", ridge_N);
	
	int i;  for(i = 0; i<vertex_N; i++)
		if(mesh->isRidge[i]) {
			fprintf(file, "%f ", vertex[vertex_table[i]][0]);
			fprintf(file, "%f ", vertex[vertex_table[i]][1]);
			fprintf(file, "%f\n", vertex[vertex_table[i]][2]);
		}
	fprintf(file, "%d\n", ravine_N);
	int i;  for(i = 0; i<vertex_N; i++)
		if(mesh->isRavine[i]) {
			fprintf(file, "%f ", vertex[vertex_table[i]][0]);
			fprintf(file, "%f ", vertex[vertex_table[i]][1]);
			fprintf(file, "%f\n", vertex[vertex_table[i]][2]);
		}
	delete[] vertex_table;
	//保存谷脊线数目及坐标结束
	osg::ref_ptr<osg::Node> meshModel = new osg::Node();
	return meshModel;
}