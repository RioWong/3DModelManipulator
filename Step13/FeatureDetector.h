
#pragma once


#include "MeshData.h"
#include "math.h"
#include <float.h>
#include "Node.h"
#include "MyApp.h"
class FeatureDetector  
{
public:
	MeshData* mesh;
	Node *ridge_pair, *ravine_pair, *normal_pair;
	int *ridge_path_N, *ravine_path_N, *normal_path_N;
	Node ***ridge_path, ***ravine_path, ***normal_path;

public:
	void ridgeBelyaevCGI98S();
	void ridgeBelyaevCGI98();
	void ridgeTriangle();
	void connectRRedge(float T, double T_ridge, double T_ravine);
	void eliminateBranchRavine(int from, int to);
	double checkBranchRavine(int index, double K, BOOL *isVisit, double T);
	void eliminateBranchRidge(int from, int to);
	double checkBranchRidge(int index, double K, BOOL* isVisit, double T);
	void connectRavineBF(int i, double K);
	void connectRidgeBF(int i, double K);
	double connectRavineDF(int i, double K);
	double connectRidgeDF(int i, double K);
	void quickSortD(int *index, double *w, int start, int end);
	void setRidgeEdges();
	void addBridgeRR(double T1, double T2);
	Node* searchPathRavine(int start, float T, double *k, float K, int *degree, int *id, double KT, BOOL isCross, double* K2);
	Node* searchPathRidge(int start, float T, double *k, float K, int *degree, int *id, double KT, BOOL isCross, double *K2);
	void connectRR(float T, double T1, double T2, BOOL isCross);
	void crossThresholding();
	void automaticThresholding(float hi, float low, double& ridge_low, double& ravine_low);
	void generateRidgeLine2(float length, float step, double T1, double T2);
	void generateRidgeLine(float length, float step, double T1, double T2);
	void subPixel();
	void ridgeT();
	void detectAngleMax(float hi, float row);
	void quickSortD(double *w, int start, int end);
	void detectShapEdge(float hi, float row);
	void setRREdge2(float hi, float low, float T);
	void setRREdge(float hi, float low, float T);
	void quickSort(int *index, float *w, int start, int end);
	void connectRRE2(float T);
	Node* searchPath2(int start, float T, float **weight, int *degree, int *id);
	void connectRRE(float T);
	Node* searchPath(int start, float T, double *k, float K, int *degree, int* id);
	void connectRR(float T);
	Node* computePath(int start, int end, double* k, float K, BOOL *valid_vertex);
	void connectPairVertex(Node** ridge_cluster, Node** ravine_cluster);
	void generatePairVertexRavine(Node** ravine_edge);
	void generatePairVertexRidge(Node **ridge_edge);
	void quickSort(float* w, int start, int end);
	void generatePairVertex(Node** ridge_edge, Node** ravine_edge);
	void generateRavinePairVertex();
	void connectPairVertex();
	Node** computePath(int start, int* ends, int end_N, double* k, float K);
	void generateRidgePairVertex();
	void precisionRidge();
	void automaticThresholding(float hi, float low);
	void setMeshData(MeshData* mesh);
	void constantThresholding(float T);
	void ridge();
	FeatureDetector();
	virtual ~FeatureDetector();

protected:
	inline void upheap(float *a, int N, int k, int *p, int *q);
	inline void downheap(float *a, int N, int k, int *p, int *q);
private:
	void labelConnectedComponent(int v, int initial, int current_id, int *id);
};
