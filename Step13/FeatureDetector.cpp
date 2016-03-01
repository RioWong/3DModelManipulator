/************************************************************************
Yutaka Ohtake 

FeatureDetector.cpp
Extraction of curvature features.

Copyright (c) 1999-2001 
The University of Aizu. All Rights Reserved.
************************************************************************/

//#include "stdafx.h"
//#include "MeshEditor.h"
#include "FeatureDetector.h"



#define PI 3.14159265

//////////////////////////////////////////////////////////////////////
// \’z/Á–Å
//////////////////////////////////////////////////////////////////////

FeatureDetector::FeatureDetector()
{
	ridge_pair = NULL;
	ravine_pair = NULL;
	normal_pair = NULL;
	mesh = NULL;
}

FeatureDetector::~FeatureDetector()
{

}

void FeatureDetector::ridge()
{
	int i;
	int vertex_N = mesh->vertex_N;
	//allocate memory
	if(mesh->isRidge != NULL)
		delete[] mesh->isRidge;
	mesh->isRidge = new BOOL[vertex_N];
	BOOL *isRidge = mesh->isRidge;

	if(mesh->isRavine != NULL)
		delete[] mesh->isRavine;
	mesh->isRavine = new BOOL[vertex_N];
	BOOL *isRavine = mesh->isRavine;

	if(mesh->k_max == NULL || mesh->k_min == NULL 
		|| mesh->t_max == NULL || mesh->t_min == NULL)
		mesh->computePrincipal();

	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	double (*t_max)[3] = mesh->t_max;
	double (*t_min)[3] = mesh->t_min;

	if(mesh->normal == NULL)
		mesh->computeNormal();

	float (*normal)[3] = mesh->normal;

	float (*vertex)[3] = mesh->vertex;
	
	for(i=0; i<vertex_N; i++){
		//skip boundary points
		if(mesh->isBound[i]){
			isRidge[i] = isRavine[i] = false;
			continue;
		}

		int nei_N, *nei;
		mesh->getConsistent1Ring(i, nei, nei_N);
		if(nei_N == 0){
			isRidge[i] = isRavine[i] = false;
			continue;
		}
	
		double k_max_max1, k_max_max2, k_max_min1, k_max_min2;
		double k_min_max1, k_min_max2, k_min_min1, k_min_min2;
		double l_max1, l_max2, l_min1, l_min2;

		double n1[3];
		MeshData::CROSS(n1, mesh->normal[i], mesh->t_max[i]);
		double n2[3];
		MeshData::CROSS(n2, mesh->normal[i], mesh->t_min[i]);

		for(int m=0; m<nei_N; m++){
			int j = nei[m];
			int k = nei[(m+1)%nei_N];

			float t1[3];
			MeshData::VEC(t1, vertex[j], vertex[i]);
			float t2[3];
			MeshData::VEC(t2, vertex[k], vertex[i]);
			
			double in1 = MeshData::DOT(n1, t1);
			double in2 = MeshData::DOT(n1, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = MeshData::AREA(o, t1, v);
				double area2 = MeshData::AREA(o, t2, v);
				double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
				double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

				if(isOne){
					k_max_max1 = k_max_in;
					k_min_max1 = k_min_in;
					l_max1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
				else{
					k_max_max2 = k_max_in;
					k_min_max2 = k_min_in;
					l_max2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
			}

			in1 = MeshData::DOT(n2, t1);
			in2 = MeshData::DOT(n2, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = MeshData::AREA(o, t1, v);
				double area2 = MeshData::AREA(o, t2, v);
				double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
				double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

				if(isOne){
					k_max_min1 = k_max_in;
					k_min_min1 = k_min_in;
					l_min1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}
				else{
					k_max_min2 = k_max_in;
					k_min_min2 = k_min_in;
					l_min2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
				}	
			}
		}

		double d1 = sqrt(l_min1);
		double d2 = sqrt(l_min2);
		
		double b = (k_max_min2-k_max[i])/d2 + (k_max[i]-k_max_min1)/d1 
			- (k_max_min2-k_max_min1)/(d1+d2);

		double K1 = k_max_max1; // + 0.5*l_max1*b*b/(k_max[i]-k_min[i]);
		double K2 = k_max_max2; // + 0.5*l_max2*b*b/(k_max[i]-k_min[i]);

		isRidge[i] = ((k_max[i] > K1) && (k_max[i] > K2)); // || ((k_max[i] < K1) && (k_max[i] < K2));
	
		d1 = sqrt(l_max1);
		d2 = sqrt(l_max2);
	
		b = (k_min_max2-k_min[i])/d2 + (k_min[i]-k_min_max1)/d1 
			- (k_min_max2-k_min_max1)/(d1+d2);

		K1 = k_min_min1; // - 0.5*l_min1*b*b/(k_max[i]-k_min[i]);
		K2 = k_min_min2; // - 0.5*l_min2*b*b/(k_max[i]-k_min[i]);

		isRavine[i] = ((k_min[i] < K1) && (k_min[i] < K2)); // || ((k_min[i] > K1) && (k_min[i] > K2));
	}
}

void FeatureDetector::constantThresholding(float T)
{
	int vertex_N = mesh->vertex_N;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	int i;  for(i = 0; i<vertex_N; i++){
		isRidge[i] &= (k_max[i] > T);
		isRavine[i] &= (k_min[i] < -T);
	}
}

void FeatureDetector::setMeshData(MeshData *mesh)
{
	this->mesh = mesh;
}

void FeatureDetector::automaticThresholding(float hi, float low)
{
	int vertex_N = mesh->vertex_N;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int *isBound = mesh->isBound;

	BOOL* weak_ridge = new BOOL[vertex_N];
	BOOL* weak_ravine = new BOOL[vertex_N];
	int i;  for(i = 0; i<vertex_N; i++){
		weak_ridge[i] = false;
		weak_ravine[i] = false;
	}
/*
	float k_max_ave = 0;
	float k_min_ave = 0;

	float k_max_var = 0;
	float k_min_var = 0;

	for(i=0; i<vertex_N; i++){
		if(mesh->vertex_link[i]->next != NULL && !isBound[i]){
			k_max_ave += k_max[i];
			k_min_ave += k_min[i];
			k_max_var += k_max[i]*k_max[i];
			k_min_var += k_min[i]*k_min[i];
		}
	}
	k_max_ave /= (float)vertex_N;
	k_min_ave /= (float)vertex_N;

	k_max_var = k_max_var/((float)vertex_N) - k_max_ave*k_max_ave;
	k_min_var = k_min_var/((float)vertex_N) - k_min_ave*k_min_ave;
	
	float ridge_hi = k_max_ave + 0.25f*(float)sqrt(k_max_var);
	float ridge_ro = k_max_ave - 0.53f*(float)sqrt(k_max_var);

	float ravine_hi = k_min_ave - 0.25f*(float)sqrt(k_min_var);
	float ravine_ro = k_min_ave + 0.53f*(float)sqrt(k_min_var);
*/
	float *k1 = new float[vertex_N];
	float *k2 = new float[vertex_N];

	int c_max = 0;
	int c_min = 0;
	for(i=0; i<vertex_N; i++){
		if(!_isnan(k_max[i]) && k_max[i] > 0)
			k1[c_max++] = (float)k_max[i];
		if(!_isnan(k_min[i]) && k_min[i] < 0)
			k2[c_min++] = (float)k_min[i];
	}
	quickSort(k1, 0, c_max-1);
	quickSort(k2, 0, c_min-1);

	float ridge_hi = k1[(int)(hi*(c_max-1))];
	float ridge_ro = k1[(int)(low*(c_max-1))];

	float ravine_hi = k2[(int)((1.0-hi)*(c_min-1))];
	float ravine_ro = k2[(int)((1.0-low)*(c_min-1))];

	delete[] k1;
	delete[] k2;

	for(i=0; i<vertex_N; i++){
		if(isBound[i])
			continue;

		if(isRidge[i]){
			if(k_max[i] > ridge_hi)
				isRidge[i] = true;
			else
				isRidge[i] = false;
			if(k_max[i] > ridge_ro)
				weak_ridge[i] = true;
			else
				weak_ridge[i] = false;
		}
		else
			weak_ridge[i] = false;

		if(isRavine[i]){
			if(k_min[i] < ravine_hi)
				isRavine[i] = true;
			else
				isRavine[i] = false;
			if(k_min[i] < ravine_ro)
				weak_ravine[i] = true;
			else
				weak_ravine[i] = false;
		}
		else
			weak_ravine[i] = false;
	}

	BOOL growth_flag = true;
	while(growth_flag){
		growth_flag = false;
		for(i=0; i<vertex_N; i++){
			if(isBound[i])
				continue;

			if(isRidge[i]){
				int *nei, nei_N;
				mesh->get1Ring(i, nei, nei_N);
				
				int j;  for(i = 0; j<nei_N; j++){
					int k = nei[j];
					if(weak_ridge[k] && !isRidge[k]){
						isRidge[k] = true;
						growth_flag = true;
					}
				}
			}

			if(isRavine[i]){
				int *nei, nei_N;
				mesh->get1Ring(i, nei, nei_N);
				
				int j;  for(i = 0; j<nei_N; j++){
					int k = nei[j];
					if(weak_ravine[k] && !isRavine[k]){
						isRavine[k] = true;
						growth_flag = true;
					}
				}
			}
		}
	}
	delete[] weak_ridge;
	delete[] weak_ravine;
}

void FeatureDetector::precisionRidge()
{
	/*
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;
	float (*vertex)[3] = mesh->vertex;
	float (*normal)[3] = mesh->normal;
	float *k_max = mesh->k_max;
	float *k_min = mesh->k_min;
	float (*t_max)[3] = mesh->t_max;
	float (*t_min)[3] = mesh->t_min;
	int* isBound = mesh->isBound;
	BOOL* isRidge = mesh->isRidge;
	BOOL* isRavine = mesh->isRavine;
	int i;

	if(mesh->ridge_tri == NULL)
		mesh->ridge_tri = new BOOL[mesh->face_N];
	BOOL* ridge_tri = mesh->ridge_tri;

	if(mesh->ravine_tri == NULL)
		mesh->ravine_tri = new BOOL[mesh->face_N];
	BOOL* ravine_tri = mesh->ravine_tri;

	for(i=0; i<face_N; i++){
		ridge_tri[i] = false;
		ravine_tri[i] = false;
	}

	for(i=0; i<vertex_N; i++){
		if(isBound[i])
			continue;
		
		if(isRidge[i]){
			int nei_N, *nei;
			mesh->getConsistent1Ring(i, nei, nei_N);
			
			int faceA = -1;
			int faceB = -1;
			float PA, PB, ka, kb;
			float n[3];
			int cA = 0;
			int cB = 0;
			MeshData::CROSS(n, normal[i], t_max[i]);
			for(int m=0; m<nei_N; m++){
				int j = nei[m];
				int k = nei[(m+1)%nei_N];

				float t1[3];
				MeshData::VEC(t1, vertex[j], vertex[i]);
				float t2[3];
				MeshData::VEC(t2, vertex[k], vertex[i]);

				float in1 = MeshData::DOT(n, t1);
				float in2 = MeshData::DOT(n, t2);
				if(in1 * in2 <= 0){
					float v[3];
					BOOL isA = true;
					if(in1 < 0){
						in1 = -in1;
					}
					if(in2 < 0){
						in2 = -in2;
						isA = false;
					}
					v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
					v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
					v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

					float o[3] = {0,0,0};
					float area1 = MeshData::AREA(o, t1, v);
					float area2 = MeshData::AREA(o, t2, v);
					float k_jk = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
					
					if(isA){
						cA++;
						faceA = face_id;
						PA = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
						ka = k_max[i] - k_jk;
					}
					else{
						cB++;
						faceB = face_id;
						PB = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
						kb = k_max[i] - k_jk;
					}
				}
			}
			if(faceA >= 0 && faceB >= 0){
				if(PB*ka < PA*kb)
					ridge_tri[faceA] = true;
				else
					ridge_tri[faceB] = true;
			}
			if(faceA >= 0 && faceB < 0)
				ridge_tri[faceA] = true;
			else if(faceA < 0 && faceB >= 0)
				ridge_tri[faceB] = true;
		}

		if(isRavine[i]){
			int nei_N, *nei;
			mesh->getConsistent1Ring(i, nei, nei_N);
	
			int faceA = -1;
			int faceB = -1;
			float PA, PB, ka, kb;
			float n[3];
			MeshData::CROSS(n, normal[i], t_min[i]);
			for(int m=0; m<nei_N; m++){
				int j = nei[m];
				int k = nei[(m+1)%nei_N];

				float t1[3];
				MeshData::VEC(t1, vertex[j], vertex[i]);
				float t2[3];
				MeshData::VEC(t2, vertex[k], vertex[i]);

				float in1 = MeshData::DOT(n, t1);
				float in2 = MeshData::DOT(n, t2);
				if(in1 * in2 <= 0){
					float v[3];
					BOOL isA = true;
					if(in1 < 0){
						in1 = -in1;
					}
					if(in2 < 0){
						in2 = -in2;
						isA = false;
					}
					v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
					v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
					v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

					float o[3] = {0,0,0};
					float area1 = MeshData::AREA(o, t1, v);
					float area2 = MeshData::AREA(o, t2, v);
					float k_jk = (area2*k_min[j] + area1*k_min[k])/(area1+area2);
					
					if(isA){
						faceA = face_id;
						PA = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
						ka = k_min[i] - k_jk;
					}
					else{
						faceB = face_id;
						PB = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
						kb = k_min[i] - k_jk;
					}
				}
			}
			if(faceA >= 0 && faceB >= 0){
				if(PB*ka < PA*kb)
					ravine_tri[faceA] = true;
				else
					ravine_tri[faceB] = true;
			}
			if(faceA >= 0 && faceB < 0)
				ravine_tri[faceA] = true;
			else if(faceA < 0 && faceB >= 0)
				ravine_tri[faceB] = true;
		}
	}*/
}

void FeatureDetector::generateRidgePairVertex()
{
	int vertex_N = mesh->vertex_N;
	BOOL* isRidge = mesh->isRidge;

	if(ridge_pair != NULL)
		delete ridge_pair;

	ridge_pair = new Node;

	int i;  for(i = 0; i<vertex_N; i++){
		if(!isRidge[i])
			continue;

		int *nei, nei_N;
		mesh->getConsistent1Ring(i, nei, nei_N);
		if(nei_N == 0)
			continue;
		
		if(isRidge[i]){
			int j;  for(i = 0; j<nei_N; j++){
				if(isRidge[nei[j]]){
					ridge_pair->append(i, nei[j]);
				}
			}
		}
	}
}

Node** FeatureDetector::computePath(int start, int* ends, int end_N, double* k, float K)
{
	int vertex_N = mesh->vertex_N;
	float (*vertex)[3] = mesh->vertex;
	float* label = new float[vertex_N];
	int* pred = new int[vertex_N];

	int* heap = new int[vertex_N];
	int* index = new int[vertex_N];
	int last_heap = 0;
	int i;  for(i = 0; i<vertex_N; i++)
		index[i] = -1;

	BOOL* is_end = new BOOL[end_N];
	for(i=0; i<end_N; i++)
		is_end[i] = false;

	for(i=0; i<vertex_N; i++){
		label[i] = 1000000;
		pred[i] = i;
	}
	label[start] = 0;
	
	int u = start;
	int last_u = -1;
	while(u != last_u){
		int *nei, nei_N;
		mesh->getConsistent1Ring(u, nei, nei_N);
		for(i=0; i<nei_N; i++){
			int v = nei[i];
			float w_k = (float)(0.5*((k[u]-K)*(k[u]-K) + (k[v]-K)*(k[v]-K)) 
				+ (k[u]-k[v])*(k[u]-k[v])/6.0 + 1.0);
			float w = (float)(w_k*MeshData::DIST(vertex[u], vertex[v]));
			/*
			float abc = MeshData::DIST(vertex[u], vertex[v])
							*MeshData::DIST(vertex[v], vertex[pred[u]])
								*MeshData::DIST(vertex[pred[u]], vertex[u]);
								*/
			/*
			float w2;
			if(abc != 0)
				w2 = 4*MeshData::AREA(vertex[u], vertex[v], vertex[pred[u]])/abc;
			else
				w2 = 0;
				*/
			
			float M = label[u] + w;
			if(index[v] < 0){
				label[v] = M;
				pred[v] = u;
				//insert;
				heap[++last_heap] = v;
				index[v] = last_heap;
				//upheap;
				upheap(label, last_heap, last_heap, heap, index);
			}
			else if(M < label[v]){
				label[v] = M;
				pred[v] = u;
				//change;
				if(index[v] != 1 && M < label[heap[index[v]/2]])
					//upheap;
					upheap(label, last_heap, index[v], heap, index);
				else
					//downheap;
					downheap(label, last_heap, index[v], heap, index);
			}
		}

		for(i=0 ;i<end_N; i++)
			if(u == ends[i])
				is_end[i] = true;
		BOOL terminate = true;
		for(i=0; i<end_N; i++)
			terminate &= is_end[i];
		if(terminate) 
			break;

		//delete;
		last_u = u;
		u = heap[1];
		if(last_heap == 0)
			break;
		heap[1] = heap[last_heap--];
		index[heap[1]] = 1;
		//down heap;
		downheap(label, last_heap, 1, heap, index);
	}
	delete[] index;
	delete[] heap;
	delete[] label;

	Node** pathes = (Node**)new Node[end_N];
	for(i=0; i<end_N; i++){
		pathes[i] = new Node;
		if(!is_end[i]){
			continue;
		}
		for(int current = ends[i]; current != start; current = pred[current])
			pathes[i]->append(current, -1);
		pathes[i]->append(start, -1);
	}
	
	delete[] pred;
	delete[] is_end;

	return pathes;
}

void FeatureDetector::connectPairVertex()
{
	int vertex_N = mesh->vertex_N;

	int ends[100];
	ridge_path = (Node***)new Node[vertex_N];
	ravine_path = (Node***)new Node[vertex_N];
	normal_path = (Node***)new Node[vertex_N];

	ridge_path_N = new int[vertex_N];
	ravine_path_N = new int[vertex_N];
	normal_path_N = new int[vertex_N];
	int i;  for(i = 0; i<vertex_N; i++){
		ridge_path_N[i] = 0;
		ravine_path_N[i] = 0;
		normal_path_N[i] = 0;
	}

	Node* current = ridge_pair;
	while(true){
		if(current->next == NULL)
			break;
		int start = current->v;
		int end_N = 0;
		for(;current->v==start; current=current->next){
			ends[end_N] = current->f;
			end_N++;
		}
		ridge_path[start] = computePath(start, ends, end_N, mesh->k_max, 10);
		ridge_path_N[start] = end_N;
	}
	
	current = ravine_pair;
	while(true){
		if(current->next == NULL)
			break;
		int start = current->v;
		int end_N = 0;
		for(;current->v==start; current=current->next){
			ends[end_N] = current->f;
			end_N++;
		}
		ravine_path[start] = computePath(start, ends, end_N, mesh->k_min, -10);
		ravine_path_N[start] = end_N;
	}

	/*
	current = normal_pair;
	float* dummy = new float[vertex_N];
	for(i=0; i<vertex_N; i++)
		dummy[i] = 0;
	while(true){
		if(current->next == NULL)
			break;
		int start = current->v;
		int end_N = 0;
		for(;current->v==start; current=current->next){
			ends[end_N] = current->f;
			end_N++;
		}
		normal_path[start] = computePath(start, ends, end_N, dummy, 10);
		normal_path_N[start] = end_N;
	}*/
}

void FeatureDetector::generateRavinePairVertex()
{
	int vertex_N = mesh->vertex_N;
	BOOL* isRavine = mesh->isRavine;

	if(ravine_pair != NULL)
		delete ravine_pair;

	ravine_pair = new Node;

	int i;  for(i = 0; i<vertex_N; i++){
		if(!isRavine[i])
			continue;

		int *nei, nei_N;
		mesh->getConsistent1Ring(i, nei, nei_N);
		if(nei_N == 0)
			continue;
		
		if(isRavine[i]){
			int j;  for(i = 0; j<nei_N; j++){
				if(isRavine[nei[j]])
					ravine_pair->append(i, nei[j]);
			}
		}
	}
}

inline void FeatureDetector::upheap(float *a, int N, int k, int *p, int *q)
{
	int v;
	v = p[k];
	while(k > 1 && a[p[k/2]] >= a[v]){
		p[k] = p[k/2]; q[p[k/2]] = k; k = k/2;
	}
	p[k] = v; q[v] = k;
}

inline void FeatureDetector::downheap(float *a, int N, int k, int *p, int *q)
{
	int j, v;
	v = p[k];
	while(k <= N/2){
		j = k+k;
		if(j < N && a[p[j]] > a[p[j+1]]) j++;
		if(a[v] <= a[p[j]]) break;
		p[k] = p[j]; q[p[j]] = k; k = j;
	}
	p[k] = v; q[v] = k;
}

void FeatureDetector::generatePairVertex(Node **ridge_edge, Node **ravine_edge)
{
	int vertex_N = mesh->vertex_N;
	
	if(ridge_pair != NULL)
		delete ridge_pair;

	ridge_pair = new Node;
	
	if(ravine_pair != NULL)
		delete ravine_pair;

	ravine_pair = new Node;

	if(normal_pair != NULL)
		delete normal_pair;

	normal_pair = new Node;

	int i;  for(i = 0; i<vertex_N; i++){
		int *nei, nei_N;
		mesh->getConsistent1Ring(i, nei, nei_N);
		if(nei_N == 0)
			continue;
		int j;  for(i = 0; j<nei_N; j++){

			BOOL flag = true;
			for(Node* current = ridge_edge[i]; current->next != NULL; current = current->next)
				if(current->v == nei[j])
					flag = false;
			for(current = ravine_edge[i]; current->next != NULL; current = current->next)
				if(current->v == nei[j])
					flag = false;

			if(flag && nei[j] < i)
				normal_pair->append(i, nei[j]);
		}

		for(Node* current = ridge_edge[i]; current->next != NULL; current = current->next)
			if(current->v < i)
				ridge_pair->append(i, current->v);

		for(current = ravine_edge[i]; current->next != NULL; current = current->next)
			if(current->v < i)
				ravine_pair->append(i, current->v);
	}
}

void FeatureDetector::quickSort(float *w, int start, int end)
{
	if(start < end){
		float v = w[end];
		int i = start-1;
		int j = end;
		while(j > i){
			for(i = i+1; w[i] < v && i < end; i++);
			for(j = j-1; w[j] > v && j >= start; j--);
			float t = w[i]; 
			w[i] = w[j]; 
			w[j] = t;
		}
		float t = w[j];
		w[j] = w[i];
		w[i] = w[end];
		w[end] = t;

		quickSort(w, start, i-1);
		quickSort(w, i+1, end);
	}
	else
		return;
}

void FeatureDetector::generatePairVertexRidge(Node **ridge_edge)
{
	int vertex_N = mesh->vertex_N;
	
	if(ridge_pair != NULL)
		delete ridge_pair;

	ridge_pair = new Node;
	
	int i;  for(i = 0; i<vertex_N; i++){
		for(Node* current = ridge_edge[i]; current->next != NULL; current = current->next)
			if(current->v < i)
				ridge_pair->append(i, current->v);
	}
}

void FeatureDetector::generatePairVertexRavine(Node **ravine_edge)
{
	int vertex_N = mesh->vertex_N;
	
	if(ravine_pair != NULL)
		delete ravine_pair;

	ravine_pair = new Node;

	int i;  for(i = 0; i<vertex_N; i++){
		for(Node* current = ravine_edge[i]; current->next != NULL; current = current->next)
			if(current->v < i)
				ravine_pair->append(i, current->v);
	}
}

void FeatureDetector::connectPairVertex(Node **ridge_cluster, Node **ravine_cluster)
{
	int vertex_N = mesh->vertex_N;

	ridge_path = (Node***)new Node[vertex_N];
	ravine_path = (Node***)new Node[vertex_N];
	normal_path = (Node***)new Node[vertex_N];

	ridge_path_N = new int[vertex_N];
	ravine_path_N = new int[vertex_N];
	normal_path_N = new int[vertex_N];
	int i;  for(i = 0; i<vertex_N; i++){
		ridge_path_N[i] = 0;
		ravine_path_N[i] = 0;
		normal_path_N[i] = 0;
	}

	BOOL *valid = new BOOL[vertex_N];

	for(Node* current=ridge_pair; current->next!=NULL; current=current->next)
		ridge_path_N[current->v]++;
	for(i=0; i<vertex_N; i++)
		if(ridge_path_N[i] > 0){
			ridge_path[i] = (Node**)new Node[ridge_path_N[i]];
			ridge_path_N[i] = 0;
		}
	for(current=ridge_pair; current->next!=NULL; current=current->next){
		int start = current->v;
		int end = current->f;

		for(i=0; i<vertex_N; i++)
			valid[i] = false;
		for(Node* c = ridge_cluster[start]; c->next!=NULL; c=c->next)
			valid[c->v] = true;
		for(c = ridge_cluster[end]; c->next!=NULL; c=c->next)
			valid[c->v] = true;
		float K = -100000;
		for(i=0; i<vertex_N; i++)
			if(valid[i] && K < mesh->k_max[i])
				K = (float)mesh->k_max[i];

		ridge_path[start][ridge_path_N[start]] 
			= computePath(start, end, mesh->k_max, K, valid);

		ridge_path_N[start]++;
	}


	for(current=ravine_pair; current->next!=NULL; current=current->next)
		ravine_path_N[current->v]++;
	for(i=0; i<vertex_N; i++)
		if(ravine_path_N[i] > 0){
			ravine_path[i] = (Node**)new Node[ravine_path_N[i]];
			ravine_path_N[i] = 0;
		}
	for(current=ravine_pair; current->next!=NULL; current=current->next){
		if(current->next == NULL)
			break;
		int start = current->v;
		int end = current->f;

		for(i=0; i<vertex_N; i++)
			valid[i] = false;
		for(Node* c = ravine_cluster[start]; c->next!=NULL; c=c->next)
			valid[c->v] = true;
		for(c = ravine_cluster[end]; c->next!=NULL; c=c->next)
			valid[c->v] = true;
		float K = 100000;
		for(i=0; i<vertex_N; i++)
			if(valid[i] && K > mesh->k_min[i])
				K = (float)mesh->k_min[i];

		ravine_path[start][ravine_path_N[start]] 
			= computePath(start, end, mesh->k_min, K, valid);

		ravine_path_N[start]++;
	}
	
	delete valid;
}

Node* FeatureDetector::computePath(int start, int end, double *k, float K, BOOL *valid_vertex)
{
	int vertex_N = mesh->vertex_N;
	float (*vertex)[3] = mesh->vertex;
	float* label = new float[vertex_N];
	int* pred = new int[vertex_N];

	int* heap = new int[vertex_N];
	int* index = new int[vertex_N];
	int last_heap = 0;
	int i;  for(i = 0; i<vertex_N; i++)
		index[i] = -1;

	for(i=0; i<vertex_N; i++){
		label[i] = 1000000;
		pred[i] = i;
	}
	label[start] = 0;
	
	int u = start;
	int last_u = -1;
	while(u != last_u){
		int *nei, nei_N;
		mesh->getConsistent1Ring(u, nei, nei_N);
		for(i=0; i<nei_N; i++){
			int v = nei[i];
			if(!valid_vertex[v])
				continue;
			float w_k = (float)(0.5*((k[u]-K)*(k[u]-K) + (k[v]-K)*(k[v]-K)) 
				+ (k[u]-k[v])*(k[u]-k[v])/6.0 + 1.0);
			float w = (float)(w_k*MeshData::DIST(vertex[u], vertex[v]));
			/*
			float abc = MeshData::DIST(vertex[u], vertex[v])
							*MeshData::DIST(vertex[v], vertex[pred[u]])
								*MeshData::DIST(vertex[pred[u]], vertex[u]);
			float w2;
			if(abc != 0)
				w2 = 4*MeshData::AREA(vertex[u], vertex[v], vertex[pred[u]])/abc;
			else
				w2 = 0;*/
			
			float M = label[u] + w;
			if(index[v] < 0){
				label[v] = M;
				pred[v] = u;
				//insert;
				heap[++last_heap] = v;
				index[v] = last_heap;
				//upheap;
				upheap(label, last_heap, last_heap, heap, index);
			}
			else if(M < label[v]){
				label[v] = M;
				pred[v] = u;
				//change;
				if(index[v] != 1 && M < label[heap[index[v]/2]])
					//upheap;
					upheap(label, last_heap, index[v], heap, index);
				else
					//downheap;
					downheap(label, last_heap, index[v], heap, index);
			}
		}

		if(u == end)
			break;		

		//delete;
		last_u = u;
		u = heap[1];
		if(last_heap == 0)
			break;
		heap[1] = heap[last_heap--];
		index[heap[1]] = 1;
		//down heap;
		downheap(label, last_heap, 1, heap, index);
	}
	delete[] index;
	delete[] heap;
	delete[] label;

	Node* path = new Node;
	
	if(pred[end] == end)
		return path;

	for(int current = end; current != start; current = pred[current])
		path->append(current, -1);
	path->append(start, -1);
	
	delete[] pred;

	return path;
}

void FeatureDetector::connectRR(float T)
{
	int vertex_N = mesh->vertex_N;
	int *degree = new int[vertex_N];
	int *id = new int[vertex_N];
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;

	int i;  for(i = 0; i<vertex_N; i++){
		if(isRidge[i]){
			degree[i] = 0;
			id[i] = 0;
			int *nei, nei_N;
			mesh->getConsistent1Ring(i, nei, nei_N);
			if(nei_N == 0)
				continue;
			int j;  for(i = 0; j<nei_N; j++)
				if(isRidge[nei[j]])
					degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}

	int current_id = 1;
	for(i=0; i<vertex_N; i++)
		if(id[i] == 0){
			labelConnectedComponent(i, 0, current_id, id);
			current_id++;
		}
			

	BOOL update = true;
	while(update){
		update = false;
		for(i=0; i<vertex_N; i++){
			if(degree[i] != 0 && degree[i] != 1)
				continue;
			Node* path = searchPath(i, T, mesh->k_max, 5, degree, id);
			if(path != NULL){
				int new_id = id[path->v];
				degree[path->v]++;
				for(Node* current = path->next; current->v!=i; current=current->next){
					isRidge[current->v] = true;
					degree[current->v] = 2;
					id[current->v] = new_id;
				}
				degree[i]++;
				labelConnectedComponent(i, id[i], new_id, id);
				update = true;
				delete path;
			}
		}
		update = false;
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRidge[i] = false;

	//////////////////Ravine///////////////////
	for(i=0; i<vertex_N; i++){
		if(isRavine[i]){
			degree[i] = 0;
			id[i] = 0;
			int *nei, nei_N;
			mesh->getConsistent1Ring(i, nei, nei_N);
			if(nei_N == 0)
				continue;
			int j;  for(i = 0; j<nei_N; j++)
				if(isRavine[nei[j]])
					degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}

	current_id = 1;
	for(i=0; i<vertex_N; i++)
		if(id[i] == 0){
			labelConnectedComponent(i, 0, current_id, id);
			current_id++;
		}
			

	update = true;
	while(update){
		update = false;
		for(i=0; i<vertex_N; i++){
			if(degree[i] != 0 && degree[i] != 1)
				continue;
			Node* path = searchPath(i, T, mesh->k_min, -5, degree, id);
			if(path != NULL){
				int new_id = id[path->v];
				degree[path->v]++;
				for(Node* current = path->next; current->v!=i; current=current->next){
					isRavine[current->v] = true;
					degree[current->v] = 2;
					id[current->v] = new_id;
				}
				degree[i]++;
				labelConnectedComponent(i, id[i], new_id, id);
				update = true;
				delete path;
			}
		}
		update = false;
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRavine[i] = false;

	delete[] id;
	delete[] degree;
}

void FeatureDetector::labelConnectedComponent(int v, int initial, int current_id, int* id)
{
	id[v] = current_id;
	int *nei, nei_N;
	mesh->getConsistent1Ring(v, nei, nei_N);
	if(nei_N == 0)
		return;
	int i;  for(i = 0; i<nei_N; i++){
		if(id[nei[i]] == initial){
			labelConnectedComponent(nei[i], initial, current_id, id);
		}
	}
}

Node* FeatureDetector::searchPath(int start, float T, double *k, float K, int *degree, int *id)
{
	int vertex_N = mesh->vertex_N;
	float (*vertex)[3] = mesh->vertex;
	float* label = new float[vertex_N];
	int* pred = new int[vertex_N];

	int* heap = new int[vertex_N];
	int* index = new int[vertex_N];
	int last_heap = 0;
	int i;  for(i = 0; i<vertex_N; i++)
		index[i] = -1;

	for(i=0; i<vertex_N; i++){
		label[i] = 1000000;
		pred[i] = i;
	}
	label[start] = 0;
	
	int u = start;
	int last_u = -1;
	int end = -1;
	while(u != last_u){
		int *nei, nei_N;
		mesh->getConsistent1Ring(u, nei, nei_N);
		for(i=0; i<nei_N; i++){
			int v = nei[i];
			if(id[v] == id[start])
				continue;

			float w_k = (float)(0.5*((k[u]-K)*(k[u]-K) + (k[v]-K)*(k[v]-K)) 
				+ (k[u]-k[v])*(k[u]-k[v])/6.0 + 1.0);
			float w = (float)(w_k*MeshData::DIST(vertex[u], vertex[v]));
			
			float M = label[u] + w;
			if(index[v] < 0){
				label[v] = M;
				pred[v] = u;
				//insert;
				heap[++last_heap] = v;
				index[v] = last_heap;
				//upheap;
				upheap(label, last_heap, last_heap, heap, index);
			}
			else if(M < label[v]){
				label[v] = M;
				pred[v] = u;
				//change;
				if(index[v] != 1 && M < label[heap[index[v]/2]])
					//upheap;
					upheap(label, last_heap, index[v], heap, index);
				else
					//downheap;
					downheap(label, last_heap, index[v], heap, index);
			}
		}

		if(u == end || last_heap == 0)
			break;		

		//delete;
		last_u = u;
		u = heap[1];
		if(/*(degree[u] == 0 || degree[u] == 1)*/ degree[u] >= 0 && id[u] != id[start]){
			end = u;
			break;
		}
		else if(label[u] > T)
			break;
		else if(last_heap == 0)
			break;
		heap[1] = heap[last_heap--];
		index[heap[1]] = 1;
		//down heap;
		downheap(label, last_heap, 1, heap, index);
	}
	delete[] index;
	delete[] heap;
	delete[] label;

	
	Node* path = NULL;
	if(end != -1){
		path = new Node;
		for(int current = end; current != start; current = pred[current])
			path->append(current, -1);
		path->append(start, -1);
	}

	delete[] pred;

	return path;
}

void FeatureDetector::connectRRE(float T)
{
	int vertex_N = mesh->vertex_N;
	int *degree = new int[vertex_N];
	int *id = new int[vertex_N];
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	Node **ridge_edge = mesh->ridge_edge;
	Node **ravine_edge = mesh->ravine_edge;

	int i;  for(i = 0; i<vertex_N; i++){
		if(isRidge[i]){
			degree[i] = 0;
			id[i] = 0;
			for(Node* current=ridge_edge[i]; current->next!=NULL; current=current->next)
				degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}

	int current_id = 1;
	for(i=0; i<vertex_N; i++)
		if(id[i] == 0){
			labelConnectedComponent(i, 0, current_id, id);
			current_id++;
		}

	for(i=0; i<vertex_N; i++){
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPath(i, T, mesh->k_max, 5, degree, id);
		if(path != NULL){
			int new_id = id[path->v];
			degree[path->v]++;
			ridge_edge[path->v]->append(path->next->v, -1);
			ridge_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRidge[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ridge_edge[current->v]->append(current->next->v, -1);
				ridge_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRidge[i] = false;



	//////////////////Ravine///////////////////
	for(i=0; i<vertex_N; i++){
		if(isRavine[i]){
			degree[i] = 0;
			id[i] = 0;
			for(Node* current=ravine_edge[i]; current->next!=NULL; current=current->next)
				degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}

	current_id = 1;
	for(i=0; i<vertex_N; i++)
		if(id[i] == 0){
			labelConnectedComponent(i, 0, current_id, id);
			current_id++;
		}
			
	for(i=0; i<vertex_N; i++){
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPath(i, T, mesh->k_min, -5, degree, id);
		if(path != NULL){
			int new_id = id[path->v];
			degree[path->v]++;
			ravine_edge[path->v]->append(path->next->v, -1);
			ravine_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRavine[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ravine_edge[current->v]->append(current->next->v, -1);
				ravine_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
	}
	
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRavine[i] = false;

	delete[] id;
	delete[] degree;
}

Node* FeatureDetector::searchPath2(int start, float T, float **weight, int *degree, int *id)
{
	int vertex_N = mesh->vertex_N;
	float (*vertex)[3] = mesh->vertex;
	float* label = new float[vertex_N];
	int* pred = new int[vertex_N];

	int* heap = new int[vertex_N];
	int* index = new int[vertex_N];
	int last_heap = 0;
	int i;  for(i = 0; i<vertex_N; i++)
		index[i] = -1;

	for(i=0; i<vertex_N; i++){
		label[i] = 1000000;
		pred[i] = i;
	}
	label[start] = 0;
	
	int u = start;
	int last_u = -1;
	int end = -1;
	while(u != last_u){
		int *nei, nei_N;
		mesh->getConsistent1Ring(u, nei, nei_N);
		for(i=0; i<nei_N; i++){
			int v = nei[i];
			float w = weight[u][i]; 
			float M = label[u] + w;
			if(index[v] < 0){
				label[v] = M;
				pred[v] = u;
				//insert;
				heap[++last_heap] = v;
				index[v] = last_heap;
				//upheap;
				upheap(label, last_heap, last_heap, heap, index);
			}
			else if(M < label[v]){
				label[v] = M;
				pred[v] = u;
				//change;
				if(index[v] != 1 && M < label[heap[index[v]/2]])
					//upheap;
					upheap(label, last_heap, index[v], heap, index);
				else
					//downheap;
					downheap(label, last_heap, index[v], heap, index);
			}
		}

		if(u == end)
			break;		

		//delete;
		last_u = u;
		u = heap[1];
		if((degree[u] == 0 || degree[u] == 1) && id[u] != id[start]){
			end = u;
			break;
		}
		else if(label[u] > T)
			break;
		else if(last_heap == 0)
			break;
		heap[1] = heap[last_heap--];
		index[heap[1]] = 1;
		//down heap;
		downheap(label, last_heap, 1, heap, index);
	}
	delete[] index;
	delete[] heap;
	delete[] label;

	
	Node* path = NULL;
	if(end != -1){
		path = new Node;
		for(int current = end; current != start; current = pred[current])
			path->append(current, -1);
		path->append(start, -1);
	}

	delete[] pred;

	return path;
}

void FeatureDetector::connectRRE2(float T)
{
	int vertex_N = mesh->vertex_N;
	int *degree = new int[vertex_N];
	int *id = new int[vertex_N];
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	Node **ridge_edge = mesh->ridge_edge;
	Node **ravine_edge = mesh->ravine_edge;
	float *strength;
	int *index = new int[vertex_N];
	float *end_strength = new float[vertex_N];
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	float (*vertex)[3] = mesh->vertex;

	int i;  for(i = 0; i<vertex_N; i++){
		if(ridge_edge[i]->next == NULL)
			isRidge[i] = false;
		if(ravine_edge[i]->next == NULL)
			isRavine[i] = false;
	}

	for(i=0; i<vertex_N; i++){
		if(isRidge[i]){
			degree[i] = 0;
			id[i] = -1;
			for(Node* current=ridge_edge[i]; current->next!=NULL; current=current->next)
				degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -2;
		}
	}

	int current_id = 0;
	for(i=0; i<vertex_N; i++)
		if(id[i] == -1){
			labelConnectedComponent(i, -1, current_id, id);
			current_id++;
		}

	strength = new float[current_id];
	for(i=0; i<current_id; i++)
		strength[i] = 0;
	for(i=0; i<vertex_N; i++){
		if(id[i] < 0)
			continue;
		if(degree[i] == 0){
			strength[id[i]] = (float)(1.0/fabs(k_max[i]));
		}
		else{
			for(Node* current=ridge_edge[i]; current->next!=NULL; current=current->next)
				strength[id[i]] -= (float)(0.5*fabs(k_max[i])*MeshData::DIST(vertex[i], vertex[current->v]));
		}
	}
	int end_N = 0;
	for(i=0; i<vertex_N; i++){
		if(degree[i] == 0 || degree[i] == 1){
			end_strength[end_N] = strength[id[i]];
			index[end_N] = i;
			end_N++;
		}
	}
	quickSort(index, end_strength, 0, end_N-1);

	int j;  for(i = 0; j<end_N; j++){
		i = index[j];
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPath(i, T, mesh->k_max, 30, degree, id);
		if(path != NULL){
			int new_id = id[path->v];
			degree[path->v]++;
			ridge_edge[path->v]->append(path->next->v, -1);
			ridge_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRidge[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ridge_edge[current->v]->append(current->next->v, -1);
				ridge_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRidge[i] = false;

	delete[] strength;

	//////////////////Ravine///////////////////
	for(i=0; i<vertex_N; i++){
		if(isRavine[i]){
			degree[i] = 0;
			id[i] = -1;
			for(Node* current=ravine_edge[i]; current->next!=NULL; current=current->next)
				degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -2;
		}
	}

	current_id = 0;
	for(i=0; i<vertex_N; i++)
		if(id[i] == -1){
			labelConnectedComponent(i, -1, current_id, id);
			current_id++;
		}

	strength = new float[current_id];
	for(i=0; i<current_id; i++)
		strength[i] = 0;
	for(i=0; i<vertex_N; i++){
		if(id[i] < 0)
			continue;
		if(degree[i] == 0){
			strength[id[i]] = (float)(1.0/fabs(mesh->k_min[i]));
		}
		else{
			for(Node* current=ravine_edge[i]; current->next!=NULL; current=current->next)
				strength[id[i]] -= (float)(0.5*fabs(mesh->k_min[i])*MeshData::DIST(vertex[i], vertex[current->v]));
		}
	}
	end_N = 0;
	for(i=0; i<vertex_N; i++){
		if(degree[i] == 0 || degree[i] == 1){
			end_strength[end_N] = strength[id[i]];
			index[end_N] = i;
			end_N++;
		}
	}
	quickSort(index, end_strength, 0, end_N-1);
			
	for(j=0; j<end_N; j++){
		i = index[j];
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPath(i, T, mesh->k_min, -30, degree, id);
		if(path != NULL){
			int new_id = id[path->v];
			degree[path->v]++;
			ravine_edge[path->v]->append(path->next->v, -1);
			ravine_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRavine[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ravine_edge[current->v]->append(current->next->v, -1);
				ravine_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
	}
	
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRavine[i] = false;

	delete[] id;
	delete[] degree;
	delete[] index;
	delete[] strength;
	delete[] end_strength;
}

void FeatureDetector::quickSort(int *index, float *w, int start, int end)
{
	if(start < end){
		float v = w[end];
		int i = start-1;
		int j = end;
		while(j > i){
			for(i = i+1; w[i] < v; i++);
			for(j = j-1; w[j] > v; j--);
			float t = w[i]; 
			w[i] = w[j]; 
			w[j] = t;

			int tmp = index[i];
			index[i] = index[j];
			index[j] = tmp;
		}
		float t = w[j];
		w[j] = w[i];
		w[i] = w[end];
		w[end] = t;

		int tmp = index[j];
		index[j] = index[i];
		index[i] = index[end];
		index[end] = tmp;

		quickSort(index, w, start, i-1);
		quickSort(index, w, i+1, end);
	}
	else
		return;
}

void FeatureDetector::setRREdge(float hi, float low, float T)
{
	int vertex_N = mesh->vertex_N;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int *isBound = mesh->isBound;
	float (*vertex)[3] = mesh->vertex;

	BOOL* weak_ridge = new BOOL[vertex_N];
	BOOL* weak_ravine = new BOOL[vertex_N];
	int i;  for(i = 0; i<vertex_N; i++){
		weak_ridge[i] = false;
		weak_ravine[i] = false;
	}

	float *k1 = new float[vertex_N];
	float *k2 = new float[vertex_N];

	int c_max = 0;
	int c_min = 0;
	for(i=0; i<vertex_N; i++){
		if(k_max[i] > 0)
			k1[c_max++] = (float)k_max[i];
		if(k_min[i] < 0)
			k2[c_min++] = (float)k_min[i];
	}
	quickSort(k1, 0, c_max-1);
	quickSort(k2, 0, c_min-1);

	float ridge_hi = k1[(int)(hi*(c_max-1))];
	float ridge_ro = k1[(int)(low*(c_max-1))];

	float ravine_hi = k2[(int)((1.0-hi)*(c_min-1))];
	float ravine_ro = k2[(int)((1.0-low)*(c_max-1))];

	delete[] k1;
	delete[] k2;

	int ridge_N = 0;
	int ravine_N = 0;
	int *ridge_index = new int[vertex_N];
	int *ravine_index = new int[vertex_N];
	float *ridge_s = new float[vertex_N];
	float *ravine_s = new float[vertex_N];
	for(i=0; i<vertex_N; i++){
		if(isRidge[i]){
			if(k_max[i] > ridge_hi){
				isRidge[i] = true;
				ridge_index[ridge_N] = i;
				ridge_s[ridge_N] = (float)fabs(k_max[i]); 
				ridge_N++;
			}
			else
				isRidge[i] = false;
			if(k_max[i] > ridge_ro)
				weak_ridge[i] = true;
			else
				weak_ridge[i] = false;
		}
		else
			weak_ridge[i] = false;

		if(isRavine[i]){
			if(k_min[i] < ravine_hi){
				isRavine[i] = true;
				ravine_index[ravine_N] = i;
				ravine_s[ravine_N] = (float)fabs(k_min[i]); 
				ravine_N++;
			}
			else
				isRavine[i] = false;
			if(k_min[i] < ravine_ro)
				weak_ravine[i] = true;
			else
				weak_ravine[i] = false;
		}
		else
			weak_ravine[i] = false;
	}
	
	//allocate ridge & ravine edge 
	Node** ridge_edge = mesh->ridge_edge;
	Node** ravine_edge = mesh->ravine_edge;

	if(ridge_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ridge_edge[i];
	}
	if(ravine_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ravine_edge[i];
	}

	ridge_edge = mesh->ridge_edge = new Node*[vertex_N];
	ravine_edge = mesh->ravine_edge = new Node*[vertex_N];

	for(i=0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}

	//Ridge
	int *id = new int[vertex_N];
	int *degree = new int[vertex_N];
	quickSort(ridge_index, ridge_s, 0, ridge_N-1);
	
	int current_id = 0;
	for(i=0; i<vertex_N; i++){
		if(weak_ridge[i]){
			degree[i] = 0;
			id[i] = current_id;
			current_id++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}
	delete[] weak_ridge;

	for(int j=ridge_N-1; j>=0; j--){
		i = ridge_index[j];
		if(degree[i] != 0)// && degree[i] != 1)
			continue;
		Node* path = searchPath(i, T, k_max, 5, degree, id);
		if(path != NULL){
			int new_id = id[path->v];
			isRidge[path->v] = true;
			degree[path->v]++;
			ridge_edge[path->v]->append(path->next->v, -1);
			ridge_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRidge[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ridge_edge[current->v]->append(current->next->v, -1);
				ridge_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
		else{
			isRidge[i] = false;
			id[i] = -1;
			degree[i] = -1;
		}
	}

	BOOL *isDecided = new BOOL[vertex_N];
	for(i=0; i<vertex_N; i++){
		if(degree[i] == 0 || degree[i] == 1)
			isDecided[i] = false;
		else
			isDecided[i] = true;
	}

	BOOL isGrow = true;
	while(isGrow){
		isGrow = false;

		for(i=0; i<vertex_N; i++)
			if(id[i] >= 0)
				id[i] += vertex_N;

		current_id = 0;
		ridge_N = 0;
		for(i=0; i<vertex_N; i++){
			if(id[i] >= vertex_N && isRidge[i]){
				labelConnectedComponent(i, id[i], current_id, id);
				current_id++;
				ridge_N++;
			}
		}

		for(i=0; i<vertex_N; i++){
			if(id[i] >= vertex_N){
				labelConnectedComponent(i, id[i], current_id, id);
				current_id++;
			}
		}

		float *strength = new float[ridge_N];
		for(i=0; i<ridge_N; i++)
			strength[i] = 0;
		for(i=0; i<vertex_N; i++){
			if(isRidge[i])
				for(Node* current=ridge_edge[i]; current->next!=NULL; current=current->next)
					strength[id[i]] -= (float)(0.5*fabs(k_max[i])*MeshData::DIST(vertex[i], vertex[current->v]));
		}
	
		int end_N = 0;
		float *end_strength = new float[2*ridge_N];
		int *index = new int[2*ridge_N];
		for(i=0; i<vertex_N; i++){
			if(degree[i] == 1 && !isDecided[i]){
				end_strength[end_N] = strength[id[i]];
				index[end_N] = i;
				end_N++;
			}
		}
		delete[] strength;
		quickSort(index, end_strength, 0, end_N-1);
		delete[] end_strength;
TRACE("%d\n", end_N);
		int j;  for(i = 0; j<end_N; j++){
			i = index[j];
			if(degree[i] != 1)
				continue;
			Node* path = searchPath(i, T, mesh->k_max, 5, degree, id);
			if(path != NULL){
				int new_id = id[path->v];
				isRidge[path->v] = true;
				degree[path->v]++;
				ridge_edge[path->v]->append(path->next->v, -1);
				ridge_edge[path->next->v]->append(path->v, -1);
				for(Node* current = path->next; current->v!=i; current=current->next){
					isRidge[current->v] = true;
					degree[current->v] = 2;
					id[current->v] = new_id;
					ridge_edge[current->v]->append(current->next->v, -1);
					ridge_edge[current->next->v]->append(current->v, -1);
				}
				degree[i]++;
				labelConnectedComponent(i, id[i], new_id, id);
				delete path;

				isGrow = true;
			}
			else{
				isDecided[i] = true;
			}
		}
		delete[] index;
	}
	//ravine
	quickSort(ravine_index, ravine_s, 0, ravine_N-1);
	current_id = 0;
	for(i=0; i<vertex_N; i++){
		if(weak_ravine[i]){
			degree[i] = 0;
			id[i] = current_id;
			current_id++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}
	delete[] weak_ravine;

	for(j=ravine_N-1; j>=0; j--){
		i = ravine_index[j];
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPath(i, T, k_min, -5, degree, id);
		if(path != NULL){
			int new_id = id[path->v];
			isRavine[path->v] = true;
			degree[path->v]++;
			ravine_edge[path->v]->append(path->next->v, -1);
			ravine_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRavine[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ravine_edge[current->v]->append(current->next->v, -1);
				ravine_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
		else{
			isRavine[i] = false;
			id[i] = -1;
		}
	}

	delete[] id;
	delete[] degree;

	/*
	int vertex_N = mesh->vertex_N;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;

	Node** ridge_edge = mesh->ridge_edge;
	Node** ravine_edge = mesh->ravine_edge;

	if(ridge_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ridge_edge[i];
	}
	if(ravine_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ravine_edge[i];
	}

	ridge_edge = mesh->ridge_edge = (Node **)new Node[vertex_N];
	ravine_edge = mesh->ravine_edge = (Node **)new Node[vertex_N];

	int i;  for(i = 0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}

	int *degree = new int[vertex_N];
	for(i=0; i<vertex_N; i++)
		degree[i] = 0;
	
	int ridge_N = 0;
	int *index = new int[vertex_N];
	float *strength = new float[vertex_N];
	for(i=0; i<vertex_N; i++){
		if(isRidge[i]){
			strength[ridge_N] = mesh->k_max[i];
			index[ridge_N] = i;
			ridge_N++;
		}
	}


	for(i=0; i<vertex_N; i++){
		if(!isRidge[i] && !isRavine[i])
			continue;
		int *nei, nei_N;
		mesh->getConsistent1Ring(i, nei, nei_N);
		
		if(isRidge[i]){
			int m = 0;
			int j;  for(i = 0; j<nei_N; j++)
				if(isRidge[nei[j]])
					m++;
			if(m == 1 || m == 2)
				for(j=0; j<nei_N; j++)
					if(isRidge[nei[j]]){
						ridge_edge[i]->append(nei[j], -1);
						both_edge[i]->append(nei[j], -1);
					}
		}

		if(isRavine[i]){
			int m = 0;
			int j;  for(i = 0; j<nei_N; j++)
				if(isRavine[nei[j]])
					m++;
			if(m == 1 || m == 2)
				for(j=0; j<nei_N; j++)
					if(isRavine[nei[j]]){
						ravine_edge[i]->append(nei[j], -1);
						if(!isRidge[i] || !isRidge[nei[j]])
							both_edge[i]->append(nei[j], -1);
					}
		}
	}
	*/
}

void FeatureDetector::setRREdge2(float hi, float low, float T)
{
	int vertex_N = mesh->vertex_N;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int *isBound = mesh->isBound;
	float (*vertex)[3] = mesh->vertex;

	BOOL* weak_ridge = new BOOL[vertex_N];
	BOOL* weak_ravine = new BOOL[vertex_N];
	int i;  for(i = 0; i<vertex_N; i++){
		weak_ridge[i] = false;
		weak_ravine[i] = false;
	}

	float *k1 = new float[vertex_N];
	float *k2 = new float[vertex_N];

	int c_max = 0;
	int c_min = 0;
	for(i=0; i<vertex_N; i++){
		if(k_max[i] > 0)
			k1[c_max++] = (float)k_max[i];
		if(k_min[i] < 0)
			k2[c_min++] = (float)k_min[i];
	}
	quickSort(k1, 0, c_max-1);
	quickSort(k2, 0, c_min-1);

	float ridge_hi = k1[(int)(hi*(c_max-1))];
	float ridge_ro = k1[(int)(low*(c_max-1))];

	float ravine_hi = k2[(int)((1.0-hi)*(c_min-1))];
	float ravine_ro = k2[(int)((1.0-low)*(c_max-1))];

	delete[] k1;
	delete[] k2;

	for(i=0; i<vertex_N; i++){
		if(isRidge[i]){
			if(k_max[i] > ridge_hi){
				isRidge[i] = true;
			}
			else
				isRidge[i] = false;
			if(k_max[i] > ridge_ro)
				weak_ridge[i] = true;
			else
				weak_ridge[i] = false;
		}
		else
			weak_ridge[i] = false;

		if(isRavine[i]){
			if(k_min[i] < ravine_hi){
				isRavine[i] = true;
			}
			else
				isRavine[i] = false;
			if(k_min[i] < ravine_ro)
				weak_ravine[i] = true;
			else
				weak_ravine[i] = false;
		}
		else
			weak_ravine[i] = false;
	}

	BOOL *ridge_tmp = new BOOL[vertex_N];
	BOOL *ravine_tmp = new BOOL[vertex_N];
	for(i=0; i<0; i++){
		int j;  for(i = 0; j<vertex_N; j++){
			ridge_tmp[j] = false;
			ravine_tmp[j] = false;
		}

		for(j=0; j<vertex_N; j++){
			if(weak_ridge[j]){
				int *nei, nei_N;
				mesh->getConsistent1Ring(j, nei, nei_N);
				if(nei_N == 0)
					continue;
				for(int k=0; k<nei_N; k++)
					ridge_tmp[nei[k]] = true;
			}
			if(weak_ravine[j]){
				int *nei, nei_N;
				mesh->getConsistent1Ring(j, nei, nei_N);
				if(nei_N == 0)
					continue;
				for(int k=0; k<nei_N; k++)
					ravine_tmp[nei[k]] = true;
			}
		}

		for(j=0; j<vertex_N; j++){
			weak_ridge[j] |= ridge_tmp[j];
			weak_ravine[j] |= ravine_tmp[j];
		}
	}
	delete[] ridge_tmp;
	delete[] ravine_tmp;

	for(i=0; i<vertex_N; i++){
		if(weak_ridge[i])
			isRidge[i] = true;
		if(weak_ravine[i])
			isRavine[i] = true;
	}
	
	//allocate ridge & ravine edge 
	Node** ridge_edge = mesh->ridge_edge;
	Node** ravine_edge = mesh->ravine_edge;

	if(ridge_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ridge_edge[i];
	}
	if(ravine_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ravine_edge[i];
	}

	ridge_edge = mesh->ridge_edge = (Node **)new Node[vertex_N];
	ravine_edge = mesh->ravine_edge = (Node **)new Node[vertex_N];

	for(i=0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}
}

void FeatureDetector::detectShapEdge(float hi, float row)
{
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;
	int (*face)[3] = mesh->face;
	float (*vertex)[3] = mesh->vertex;
	float (*normal)[3] = mesh->normal_f;
	int (*adjacent)[3] = mesh->face_link_E;

	int convex_N = 0;
	int concave_N = 0;
	int i;  for(i = 0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;
		int j;  for(i = 0; j<3; j++){
			int pair = adjacent[i][j];
			if(i<pair){
				float c1[3];
				c1[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3;
				c1[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3;
				c1[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3;

				float c2[3];
				c2[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3;
				c2[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3;
				c2[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3;

				float v[3];
				MeshData::VEC(v, c1, c2);
				double len = MeshData::LENGTH(v);
				double dot1, dot2;
				if((float)len != 0){
					dot1 = MeshData::DOT(v, normal[i])/len;
					dot2 = -MeshData::DOT(v, normal[pair])/len;
				}
				else
					continue;

				if(dot1 > 1.0)
					dot1 = 1.0;
				else if(dot1 < -1.0)
					dot1 = -1.0;
				if(dot2 > 1.0)
					dot2 = 1.0;
				else if(dot2 < -1.0)
					dot2 = -1.0;
				double angle = acos(dot1) + acos(dot2);
				if(angle < 0 || angle > 2.0*PI)
					continue;
				if(angle > PI)
					convex_N++;
				else if(angle < PI)
					concave_N++;
			}
		}
	}

	float* ridge_w = new float[convex_N];
	float* ravine_w = new float[concave_N];
	//double* ridge_w = new double[convex_N];
	//double* ravine_w = new double[concave_N];
	int convex_i = 0;
	int concave_i = 0;
	for(i=0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;
		int j;  for(i = 0; j<3; j++){
			int pair = adjacent[i][j];
			if(i<pair){
				float c1[3];
				c1[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3;
				c1[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3;
				c1[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3;

				float c2[3];
				c2[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3;
				c2[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3;
				c2[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3;

				float v[3];
				MeshData::VEC(v, c1, c2);
				double len = MeshData::LENGTH(v);
				double dot1, dot2;
				if((float)len != 0){
					dot1 = MeshData::DOT(v, normal[i])/len;
					dot2 = -MeshData::DOT(v, normal[pair])/len;
				}
				else
					continue;
				if(dot1 > 1.0)
					dot1 = 1.0;
				else if(dot1 < -1.0)
					dot1 = -1.0;
				if(dot2 > 1.0)
					dot2 = 1.0;
				else if(dot2 < -1.0)
					dot2 = -1.0;
				double w = MeshData::DOT(normal[i], normal[pair]);
				if(w > 1.0)
					w = 1.0;
				else if(w < -1.0)
					w = -1.0;
				double angle = acos(dot1) + acos(dot2);
				if(angle < 0 || angle > 2.0*PI)
					continue;
				if(angle > PI){
					ridge_w[convex_i] = (float)w;
					convex_i++;
				}
				else if(angle < PI){
					ravine_w[concave_i] = (float)w;
					concave_i++;
				}
			}
		}
	}
	
	if(convex_N != 0){
		quickSort(ridge_w, 0, convex_N-1);
	}
	if(concave_N != 0){
		quickSort(ravine_w, 0, concave_N-1);
	}

	double ridge_hi = ridge_w[(int)((1.0-hi)*convex_N)];
	double ridge_row = ridge_w[(int)((1.0-row)*convex_N)];
	double ravine_hi = ravine_w[(int)((1.0-hi)*concave_N)];
	double ravine_row = ravine_w[(int)((1.0-row)*concave_N)];

	if(convex_N != 0)
		delete[] ridge_w;
	if(concave_N != 0)
		delete[] ravine_w;

	Node** ridge_edge = mesh->ridge_edge;
	Node** ravine_edge = mesh->ravine_edge;

	if(ridge_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ridge_edge[i];
		delete[] ridge_edge;
	}
	if(ravine_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ravine_edge[i];
		delete[] ravine_edge;
	}

	mesh->ridge_edge = new Node*[vertex_N];
	ridge_edge = mesh->ridge_edge;
	mesh->ravine_edge = new Node*[vertex_N];
	ravine_edge = mesh->ravine_edge;
	for(i=0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}

	Node** weak_ridge = new Node*[vertex_N];
	Node** weak_ravine = new Node*[vertex_N];
	for(i=0; i<vertex_N; i++){
		weak_ridge[i] = new Node;
		weak_ravine[i] = new Node;
	}

	for(i=0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;

		int j;  for(i = 0; j<3; j++){
			int pair = adjacent[i][j];
			if(i<pair){
				float c1[3];
				c1[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3;
				c1[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3;
				c1[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3;

				float c2[3];
				c2[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3;
				c2[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3;
				c2[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3;

				float v[3];
				MeshData::VEC(v, c1, c2);
				double len = MeshData::LENGTH(v);
				double dot1, dot2;
				if((float)len != 0){
					dot1 = MeshData::DOT(v, normal[i])/len;
					dot2 = -MeshData::DOT(v, normal[pair])/len;
				}
				else
					continue;
				double w = MeshData::DOT(normal[i], normal[pair]);
				if(dot1 > 1.0)
					dot1 = 1.0;
				else if(dot1 < -1.0)
					dot1 = -1.0;
				if(dot2 > 1.0)
					dot2 = 1.0;
				else if(dot2 < -1.0)
					dot2 = -1.0;
				if(acos(dot1) + acos(dot2) > PI){
					if(w < ridge_hi){
						ridge_edge[face[i][j]]->append(face[i][(j+1)%3], -1);
						ridge_edge[face[i][(j+1)%3]]->append(face[i][j], -1);
					}
					if(w < ridge_row){
						weak_ridge[face[i][j]]->append(face[i][(j+1)%3], -1);
						weak_ridge[face[i][(j+1)%3]]->append(face[i][j], -1);
					}
				}
				else if(acos(dot1) + acos(dot2) < PI){
					if(w < ravine_hi){
						ravine_edge[face[i][j]]->append(face[i][(j+1)%3], -1);
						ravine_edge[face[i][(j+1)%3]]->append(face[i][j], -1);
					}
					if(w < ravine_row){
						weak_ravine[face[i][j]]->append(face[i][(j+1)%3], -1);
						weak_ravine[face[i][(j+1)%3]]->append(face[i][j], -1);
					}
				}
			}
		}
	}

	bool isGrow = true;
	while(isGrow){
		isGrow = false;
		int i;  for(i = 0; i<vertex_N; i++){
			if(ridge_edge[i]->next != NULL){
				for(Node* current=weak_ridge[i]; current->next!=NULL; current=current->next){
					bool isDetect = false;
					for(Node* search=ridge_edge[i]; search->next!=NULL; search=search->next){
						if(current->v == search->v){
							isDetect = true;
							break;
						}
					}
					if(isDetect)
						continue;
					ridge_edge[i]->append(current->v, -1);
					ridge_edge[current->v]->append(i, -1);
					isGrow = true;
				}
			}
			if(ravine_edge[i]->next != NULL){
				for(Node* current=weak_ravine[i]; current->next!=NULL; current=current->next){
					bool isDetect = false;
					for(Node* search=ravine_edge[i]; search->next!=NULL; search=search->next){
						if(current->v == search->v){
							isDetect = true;
							break;
						}
					}
					if(isDetect)
						continue;
					ravine_edge[i]->append(current->v, -1);
					ravine_edge[current->v]->append(i, -1);
					isGrow = true;
				}
			}
		}
	}

	for(i=0; i<vertex_N; i++){
		delete weak_ridge[i];
		delete weak_ravine[i];
	}
	delete[] weak_ridge;
	delete[] weak_ravine;
}

void FeatureDetector::quickSortD(double *w, int start, int end)
{
	if(start < end){
		double v = w[end];
		int i = start-1;
		int j = end;
		while(j > i){
			for(i = i+1; w[i] < v; i++);
			for(j = j-1; w[j] > v; j--);
			double t = w[i]; 
			w[i] = w[j]; 
			w[j] = t;
		}
		double t = w[j];
		w[j] = w[i];
		w[i] = w[end];
		w[end] = t;

		quickSortD(w, start, i-1);
		quickSortD(w, i+1, end);
	}
	else
		return;
}

void FeatureDetector::detectAngleMax(float hi, float row)
{
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;
	int (*face)[3] = mesh->face;
	float (*vertex)[3] = mesh->vertex;
	float (*normal)[3] = mesh->normal_f;
	int (*adjacent)[3] = mesh->face_link_E;
	int *isBound = mesh->isBound;

	int convex_N = 0;
	int concave_N = 0;
	int i;  for(i = 0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;
		int j;  for(i = 0; j<3; j++){
			int pair = adjacent[i][j];
			if(i<pair){
				float c1[3];
				c1[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3;
				c1[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3;
				c1[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3;

				float c2[3];
				c2[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3;
				c2[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3;
				c2[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3;

				float v[3];
				MeshData::VEC(v, c1, c2);
				double len = MeshData::LENGTH(v);
				double dot1, dot2;
				if((float)len != 0){
					dot1 = MeshData::DOT(v, normal[i])/len;
					dot2 = -MeshData::DOT(v, normal[pair])/len;
				}
				else
					continue;

				if(dot1 > 1.0)
					dot1 = 1.0;
				else if(dot1 < -1.0)
					dot1 = -1.0;
				if(dot2 > 1.0)
					dot2 = 1.0;
				else if(dot2 < -1.0)
					dot2 = -1.0;
				double angle = acos(dot1) + acos(dot2);
				if(angle < 0 || angle > 2.0*PI)
					continue;
				if(angle > PI)
					convex_N++;
				else if(angle < PI)
					concave_N++;
			}
		}
	}

	float* ridge_w = new float[convex_N];
	float* ravine_w = new float[concave_N];
	//double* ridge_w = new double[convex_N];
	//double* ravine_w = new double[concave_N];
	int convex_i = 0;
	int concave_i = 0;
	for(i=0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;
		int j;  for(i = 0; j<3; j++){
			int pair = adjacent[i][j];
			if(i<pair){
				float c1[3];
				c1[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3;
				c1[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3;
				c1[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3;

				float c2[3];
				c2[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3;
				c2[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3;
				c2[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3;

				float v[3];
				MeshData::VEC(v, c1, c2);
				double len = MeshData::LENGTH(v);
				double dot1, dot2;
				if((float)len != 0){
					dot1 = MeshData::DOT(v, normal[i])/len;
					dot2 = -MeshData::DOT(v, normal[pair])/len;
				}
				else
					continue;
				if(dot1 > 1.0)
					dot1 = 1.0;
				else if(dot1 < -1.0)
					dot1 = -1.0;
				if(dot2 > 1.0)
					dot2 = 1.0;
				else if(dot2 < -1.0)
					dot2 = -1.0;
				double w = MeshData::DOT(normal[i], normal[pair]);
				if(w > 1.0)
					w = 1.0;
				else if(w < -1.0)
					w = -1.0;
				double angle = acos(dot1) + acos(dot2);
				if(angle < 0 || angle > 2.0*PI)
					continue;
				if(angle > PI){
					ridge_w[convex_i] = (float)w;
					convex_i++;
				}
				else if(angle < PI){
					ravine_w[concave_i] = (float)w;
					concave_i++;
				}
			}
		}
	}
	
	if(convex_N != 0){
		quickSort(ridge_w, 0, convex_N-1);
	}
	if(concave_N != 0){
		quickSort(ravine_w, 0, concave_N-1);
	}

	double ridge_hi = ridge_w[(int)((1.0-hi)*convex_N)];
	double ridge_row = ridge_w[(int)((1.0-row)*convex_N)];
	double ravine_hi = ravine_w[(int)((1.0-hi)*concave_N)];
	double ravine_row = ravine_w[(int)((1.0-row)*concave_N)];

	if(convex_N != 0)
		delete[] ridge_w;
	if(concave_N != 0)
		delete[] ravine_w;

	Node** ridge_edge = mesh->ridge_edge;
	Node** ravine_edge = mesh->ravine_edge;

	if(ridge_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ridge_edge[i];
		delete[] ridge_edge;
	}
	if(ravine_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ravine_edge[i];
		delete[] ravine_edge;
	}

	mesh->ridge_edge = new Node*[vertex_N];
	ridge_edge = mesh->ridge_edge;
	mesh->ravine_edge = new Node*[vertex_N];
	ravine_edge = mesh->ravine_edge;
	for(i=0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}

	Node** weak_ridge = new Node*[vertex_N];
	Node** weak_ravine = new Node*[vertex_N];
	for(i=0; i<vertex_N; i++){
		weak_ridge[i] = new Node;
		weak_ravine[i] = new Node;
	}

	for(i=0; i<face_N; i++){
		if(face[i][0] < 0)
			continue;

		int j;  for(i = 0; j<3; j++){
			int pair = adjacent[i][j];
			if(i > pair)
				continue;

			float c1[3];
			c1[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3;
			c1[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3;
			c1[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3;

			float c2[3];
			c2[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3;
			c2[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3;
			c2[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3;

			float v[3];
			MeshData::VEC(v, c1, c2);
			double len = MeshData::LENGTH(v);
			double dot1, dot2;
			if((float)len != 0){
				dot1 = MeshData::DOT(v, normal[i])/len;
				dot2 = -MeshData::DOT(v, normal[pair])/len;
			}
			else
				continue;

			if(dot1 > 1.0)
				dot1 = 1.0;
			else if(dot1 < -1.0)
				dot1 = -1.0;
			if(dot2 > 1.0)
				dot2 = 1.0;
			else if(dot2 < -1.0)
				dot2 = -1.0;
			double w = MeshData::DOT(normal[i], normal[pair]);
			bool convex;
			if(acos(dot1) + acos(dot2) > PI){
				if(w > ridge_row)
					continue;
				convex = true;
			}
			else if(acos(dot1) + acos(dot2) < PI){
				if(w > ravine_row)
					continue;
				convex = false;
			}
			else
				continue;

			int i1 = face[i][j];
			int i2 = face[i][(j+1)%3];
			int pair_j;
			if(face[pair][0] == i2)
				pair_j = 0;
			else if(face[pair][1] == i2)
				pair_j = 1;
			else if(face[pair][2] == i2)
				pair_j = 2;
			else{
				TRACE("pair_j is not found.\n");
			}
			int i3 = face[i][(j+2)%3];
			int i4 = face[pair][(pair_j+2)%3];
			
			int f1 = adjacent[i][(j+1)%3];
			int f2 = adjacent[i][(j+2)%3];
			int f3 = adjacent[pair][(pair_j+1)%3];
			int f4 = adjacent[pair][(pair_j+2)%3];

			double w1 = 0;
			int valid = 0;
			if(f1 >= 0){
				int pair_f1;
				if(face[f1][0] == i3)
					pair_f1 = adjacent[f1][2];
				else if(face[f1][1] == i3)
					pair_f1 = adjacent[f1][0];
				else
					pair_f1 = adjacent[f1][1];

				if(pair_f1 >= 0){
					valid++;
					w1 += MeshData::DOT(normal[f1], normal[pair_f1]);
				}
			}
			if(f2 >= 0){
				int pair_f2;
				if(face[f2][0] == i3)
					pair_f2 = adjacent[f2][0];
				else if(face[f2][1] == i3)
					pair_f2 = adjacent[f2][1];
				else
					pair_f2 = adjacent[f2][2];

				if(pair_f2 >= 0){
					valid++;
					w1 += MeshData::DOT(normal[f2], normal[pair_f2]);
				}
			}
			if(valid > 0)
				w1 /= (double)valid;
			else
				continue;

			double w2 = 0;
			valid = 0;
			if(f3 >= 0){
				int pair_f3;
				if(face[f3][0] == i4)
					pair_f3 = adjacent[f3][2];
				else if(face[f3][1] == i4)
					pair_f3 = adjacent[f3][0];
				else
					pair_f3 = adjacent[f3][1];

				if(pair_f3 >= 0){
					valid++;
					w2 += MeshData::DOT(normal[f3], normal[pair_f3]);
				}
			}
			if(f4 >= 0){
				int pair_f4;
				if(face[f4][0] == i4)
					pair_f4 = adjacent[f4][0];
				else if(face[f4][1] == i4)
					pair_f4 = adjacent[f4][1];
				else
					pair_f4 = adjacent[f4][2];

				if(pair_f4 >= 0){
					valid++;
					w2 += MeshData::DOT(normal[f4], normal[pair_f4]);
				}
			}
			if(valid > 0)
				w2 /= (double)valid;
			else
				continue;

			if(w < w1 && w < w2){
				if(convex){
					if(w < ridge_hi){
						ridge_edge[i1]->append(i2, -1);
						ridge_edge[i2]->append(i1, -1);
					}
					weak_ridge[i1]->append(i2, -1);
					weak_ridge[i2]->append(i1, -1);
				}
				else{
					if(w < ravine_hi){
						ravine_edge[i1]->append(i2, -1);
						ravine_edge[i2]->append(i1, -1);
					}
					weak_ravine[i1]->append(i2, -1);
					weak_ravine[i2]->append(i1, -1);
				}
			}
		}
	}

	bool isGrow = true;
	while(isGrow){
		isGrow = false;
		int i;  for(i = 0; i<vertex_N; i++){
			if(ridge_edge[i]->next != NULL){
				for(Node* current=weak_ridge[i]; current->next!=NULL; current=current->next){
					bool isDetect = false;
					for(Node* search=ridge_edge[i]; search->next!=NULL; search=search->next){
						if(current->v == search->v){
							isDetect = true;
							break;
						}
					}
					if(isDetect)
						continue;
					ridge_edge[i]->append(current->v, -1);
					ridge_edge[current->v]->append(i, -1);
					isGrow = true;
				}
			}
			if(ravine_edge[i]->next != NULL){
				for(Node* current=weak_ravine[i]; current->next!=NULL; current=current->next){
					bool isDetect = false;
					for(Node* search=ravine_edge[i]; search->next!=NULL; search=search->next){
						if(current->v == search->v){
							isDetect = true;
							break;
						}
					}
					if(isDetect)
						continue;
					ravine_edge[i]->append(current->v, -1);
					ravine_edge[current->v]->append(i, -1);
					isGrow = true;
				}
			}
		}
	}

	for(i=0; i<vertex_N; i++){
		delete weak_ridge[i];
		delete weak_ravine[i];
	}
	delete[] weak_ridge;
	delete[] weak_ravine;
}


void FeatureDetector::ridgeT()
{
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;
	//allocate memory
	if(mesh->isRidge_T != NULL)
		delete[] mesh->isRidge_T;
	mesh->isRidge_T = new BOOL[face_N];
	BOOL *isRidge_T = mesh->isRidge_T;

	if(mesh->ridge_T != NULL)
		delete[] mesh->ridge_T;
	mesh->ridge_T = new double[face_N][3];
	double (*ridge_T)[3] = mesh->ridge_T; 

	if(mesh->isRavine_T != NULL)
		delete[] mesh->isRavine_T;
	mesh->isRavine_T = new BOOL[face_N];
	BOOL *isRavine_T = mesh->isRavine_T;

	if(mesh->ravine_T != NULL)
		delete[] mesh->ravine_T;
	mesh->ravine_T = new double[face_N][3];
	double (*ravine_T)[3] = mesh->ravine_T; 

	if(mesh->k_max == NULL || mesh->k_min == NULL 
		|| mesh->t_max == NULL || mesh->t_min == NULL)
		mesh->computePrincipal();

	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	double (*t_max)[3] = mesh->t_max;
	double (*t_min)[3] = mesh->t_min;

	if(mesh->normal == NULL)
		mesh->computeNormal();

	float (*normal)[3] = mesh->normal;
	float (*vertex)[3] = mesh->vertex;
	int (*face)[3] = mesh->face;
	
	for(int t=0; t<face_N; t++){
		double kt_max[3], kt_min[3];
		for(int r=0; r<3; r++){
			int i = face[t][r];

			int nei_N, *nei;
			mesh->getConsistent1Ring(i, nei, nei_N);
			
			double k_max_max1, k_max_max2, k_max_min1, k_max_min2;
			double k_min_max1, k_min_max2, k_min_min1, k_min_min2;
			double l_max1, l_max2, l_min1, l_min2;

			double n1[3];
			MeshData::CROSS(n1, normal[i], t_max[i]);
			double n2[3];
			MeshData::CROSS(n2, normal[i], t_min[i]);

			//double *n1 = t_min[i];
			//double *n2 = t_max[i];

			for(int m=0; m<nei_N; m++){
				int j = nei[m];
				int k = nei[(m+1)%nei_N];

				float t1[3];
				MeshData::VEC(t1, vertex[j], vertex[i]);
				float t2[3];
				MeshData::VEC(t2, vertex[k], vertex[i]);
			
				double in1 = MeshData::DOT(n1, t1);
				double in2 = MeshData::DOT(n1, t2);
				if(in1 * in2 <= 0){
					double v[3];
					BOOL isOne = true;
					if(in1 < 0)
						in1 = -in1;
					if(in2 < 0){
						in2 = -in2;
						isOne = false;
					}
					v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
					v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
					v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

					float o[3] = {0,0,0};
					double area1 = MeshData::AREA(o, t1, v);
					double area2 = MeshData::AREA(o, t2, v);
					double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
					double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

					if(isOne){
						k_max_max1 = k_max_in;
						k_min_max1 = k_min_in;
						l_max1 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
					}
					else{
						k_max_max2 = k_max_in;
						k_min_max2 = k_min_in;
						l_max2 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
					}
				}

				in1 = MeshData::DOT(n2, t1);
				in2 = MeshData::DOT(n2, t2);
				if(in1 * in2 <= 0){
					double v[3];
					BOOL isOne = true;
					if(in1 < 0)
						in1 = -in1;
					if(in2 < 0){
						in2 = -in2;
						isOne = false;
					}
					v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
					v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
					v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

					float o[3] = {0,0,0};
					double area1 = MeshData::AREA(o, t1, v);
					double area2 = MeshData::AREA(o, t2, v);
					double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
					double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

					if(isOne){
						k_max_min1 = k_max_in;
						k_min_min1 = k_min_in;
						l_min1 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
					}
					else{
						k_max_min2 = k_max_in;
						k_min_min2 = k_min_in;
						l_min2 = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
					}	
				}
			}
			kt_max[r] = (k_max_max1 - k_max[i])/l_max1 +
						(k_max[i] - k_max_max2)/l_max2 -
						(k_max_max1 - k_max_max2)/(l_max1 + l_max2);
			kt_min[r] = (k_min_min1 - k_min[i])/l_min1 +
						(k_min[i] - k_min_min2)/l_min2 -
						(k_min_min1 - k_min_min2)/(l_min1 + l_min2);
			/*
			kt_max[r] = (l_max1-l_max2)*k_max[i]/(l_max1*l_max2) 
						+ l_max2*k_max_max1/(l_max1*(l_max1+l_max2)) 
						- l_max1*k_max_max2/(l_max2*(l_max1+l_max2));
			kt_min[r] = (l_min1-l_min2)*k_min[i]/(l_min1*l_min2) 
						+ l_min2*k_min_min1/(l_min1*(l_min1+l_min2)) 
						- l_min1*k_min_min2/(l_min2*(l_min1+l_min2));*/
		}
		if(MeshData::DOT(t_max[face[t][0]], t_max[face[t][1]]) < 0)
			kt_max[1] = -kt_max[1];
		if(MeshData::DOT(t_max[face[t][0]], t_max[face[t][2]]) < 0)
			kt_max[2] = -kt_max[2];

		if(MeshData::DOT(t_min[face[t][0]], t_min[face[t][1]]) < 0)
			kt_min[1] = -kt_min[1];
		if(MeshData::DOT(t_min[face[t][0]], t_min[face[t][2]]) < 0)
			kt_min[2] = -kt_min[2];

		isRidge_T[t] = false;
		isRavine_T[t] = false;
		int j;  for(i = 0; j<3; j++){
			ridge_T[t][j] = -1;
			/*
			if((kt_max[j] >= 0 && kt_max[(j+1)%3] < 0) 
				|| (kt_max[j] < 0 && kt_max[(j+1)%3] >= 0)){*/
			if((k_max[face[t][j]] >= 0.2 && k_max[face[t][(j+1)%3]] < 0.2) 
				|| (k_max[face[t][j]] < 0.2 && k_max[face[t][(j+1)%3]] >= 0.2)){
				isRidge_T[t] = true;
				ridge_T[t][j] = fabs(k_max[face[t][j]])/(fabs(k_max[face[t][j]]) + fabs(k_max[face[t][(j+1)%3]]));
				//ridge_T[t][j] = fabs(kt_max[j])/(fabs(kt_max[j])+fabs(kt_max[(j+1)%3]));
			}
			ravine_T[t][j] = -1;
			if((kt_min[j] >= 0 && kt_min[(j+1)%3] < 0) 
				|| (kt_min[j] < 0 && kt_min[(j+1)%3] >= 0)){
				isRavine_T[t] = true;
				ravine_T[t][j] = fabs(kt_min[j])/(fabs(kt_min[j])+fabs(kt_min[(j+1)%3]));
			}
		}
	}
}

void FeatureDetector::subPixel()
{
	float (*vertex)[3] = mesh->vertex;
	double* k_max = mesh->k_max;
	double* k_min = mesh->k_min;
	double (*t_max)[3] = mesh->t_max;
	double (*t_min)[3] = mesh->t_min;
	BOOL* isRidge = mesh->isRidge;
	BOOL* isRavine = mesh->isRavine;
	int ridge_point_N = 0;
	int ravine_point_N = 0;
	int vertex_N = mesh->vertex_N;
	int i;  for(i = 0; i < vertex_N; i++){
		if(isRidge[i])
			ridge_point_N++;
		if(isRavine[i])
			ravine_point_N++;
	}
	mesh->ridge_point_N = ridge_point_N;
	mesh->ravine_point_N = ravine_point_N;
	MeshData::tri_point* ridge_point = mesh->ridge_point = new MeshData::tri_point[ridge_point_N];
	MeshData::tri_point* ravine_point = mesh->ravine_point = new MeshData::tri_point[ravine_point_N];
	int ridge_N = 0;
	int ravine_N = 0;
	for(i=0; i<vertex_N; i++){
		if(!isRidge[i] && !isRavine[i])
			continue;

		int nei_N, *nei, *nei_f;
		mesh->getConsistent1Ring(i, nei, nei_f, nei_N);
		if(nei_N == 0){
			isRidge[i] = isRavine[i] = false;
			continue;
		}
	
		double k_max_max1, k_max_max2, k_max_min1, k_max_min2;
		double k_min_max1, k_min_max2, k_min_min1, k_min_min2;
		double l_max1, l_max2, l_min1, l_min2;
		int f_max1, f_max2, f_min1, f_min2;

		double n1[3];
		MeshData::CROSS(n1, mesh->normal[i], mesh->t_max[i]);
		double n2[3];
		MeshData::CROSS(n2, mesh->normal[i], mesh->t_min[i]);

		for(int m=0; m<nei_N; m++){
			int j = nei[m];
			int k = nei[(m+1)%nei_N];

			float t1[3];
			MeshData::VEC(t1, vertex[j], vertex[i]);
			float t2[3];
			MeshData::VEC(t2, vertex[k], vertex[i]);
			
			double in1 = MeshData::DOT(n1, t1);
			double in2 = MeshData::DOT(n1, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = MeshData::AREA(o, t1, v);
				double area2 = MeshData::AREA(o, t2, v);
				double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
				double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

				if(isOne){
					k_max_max1 = k_max_in;
					k_min_max1 = k_min_in;
					l_max1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
					f_max1 = nei_f[m];
				}
				else{
					k_max_max2 = k_max_in;
					k_min_max2 = k_min_in;
					l_max2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
					f_max2 = nei_f[m];
				}
			}

			in1 = MeshData::DOT(n2, t1);
			in2 = MeshData::DOT(n2, t2);
			if(in1 * in2 <= 0){
				double v[3];
				BOOL isOne = true;
				if(in1 < 0)
					in1 = -in1;
				if(in2 < 0){
					in2 = -in2;
					isOne = false;
				}
				v[0] = (in2*t1[0] + in1*t2[0])/(in1 + in2);
				v[1] = (in2*t1[1] + in1*t2[1])/(in1 + in2);
				v[2] = (in2*t1[2] + in1*t2[2])/(in1 + in2);

				float o[3] = {0,0,0};
				double area1 = MeshData::AREA(o, t1, v);
				double area2 = MeshData::AREA(o, t2, v);
				double k_max_in = (area2*k_max[j] + area1*k_max[k])/(area1+area2);
				double k_min_in = (area2*k_min[j] + area1*k_min[k])/(area1+area2);

				if(isOne){
					k_max_min1 = k_max_in;
					k_min_min1 = k_min_in;
					l_min1 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
					f_min1 = nei_f[m];
				}
				else{
					k_max_min2 = k_max_in;
					k_min_min2 = k_min_in;
					l_min2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
					f_min2 = nei_f[m];
				}	
			}
		}
		if(isRidge[i]){
			double d1 = sqrt(l_min1);
			double d2 = sqrt(l_min2);
			
			double b = (k_max_min2-k_max[i])/d2 + (k_max[i]-k_max_min1)/d1 
				- (k_max_min2-k_max_min1)/(d1+d2);

			double K1 = k_max_max1 + 0.5*l_max1*b*b/(k_max[i]-k_min[i]);
			double K2 = k_max_max2 + 0.5*l_max2*b*b/(k_max[i]-k_min[i]);

			d1 = sqrt(l_max1);
			d2 = sqrt(l_max2);

			double x = 0.5*(d2*d2*(-k_max[i]+K1) + d1*d1*(k_max[i]-K2))
				/(d2*(k_max[i]-K1) + d1*(k_max[i]-K2));
			if(x < 0)
				ridge_point[ridge_N].id = f_max1;
			else
				ridge_point[ridge_N].id = f_max2;
				
			float* nor = mesh->normal_f[ridge_point[ridge_N].id];
			double dot = MeshData::DOT(t_max[i], nor);
			double v[3];
			v[0] = t_max[i][0] - dot*nor[0];
			v[1] = t_max[i][1] - dot*nor[1];
			v[2] = t_max[i][2] - dot*nor[2];
			double len = MeshData::LENGTH(v);
			if((float)len == 0)
				v[0] = v[1] = v[2] = 0;
			else{
				v[0] = v[0]*x/len;
				v[1] = v[1]*x/len;
				v[2] = v[2]*x/len;
			}
			float p[3];
			p[0] = (float)(vertex[i][0] + v[0]);
			p[1] = (float)(vertex[i][1] + v[1]);
			p[2] = (float)(vertex[i][2] + v[2]);
			mesh->barycentricCoor(ridge_point[ridge_N].id, p, ridge_point[ridge_N].u, ridge_point[ridge_N].v);
			ridge_N++;
		}
		if(isRavine[i]){
			double d1 = sqrt(l_max1);
			double d2 = sqrt(l_max2);
	
			double b = (k_min_max2-k_min[i])/d2 + (k_min[i]-k_min_max1)/d1 
				- (k_min_max2-k_min_max1)/(d1+d2);

			double K1 = k_min_min1 - 0.5*l_min1*b*b/(k_max[i]-k_min[i]);
			double K2 = k_min_min2 - 0.5*l_min2*b*b/(k_max[i]-k_min[i]);

			d1 = sqrt(l_min1);
			d2 = sqrt(l_min2);

			double x = 0.5*(d2*d2*(-k_min[i]+K1) + d1*d1*(k_min[i]-K2))
				/(d2*(k_min[i]-K1) + d1*(k_min[i]-K2));
			if(x < 0)
				ravine_point[ravine_N].id = f_min1;
			else
				ravine_point[ravine_N].id = f_min2;
				
			float* nor = mesh->normal_f[ravine_point[ravine_N].id];
			double dot = MeshData::DOT(t_min[i], nor);
			double v[3];
			v[0] = t_min[i][0] - dot*nor[0];
			v[1] = t_min[i][1] - dot*nor[1];
			v[2] = t_min[i][2] - dot*nor[2];
			double len = MeshData::LENGTH(v);
			if((float)len == 0)
				v[0] = v[1] = v[2] = 0;
			else{
				v[0] = v[0]*x/len;
				v[1] = v[1]*x/len;
				v[2] = v[2]*x/len;
			}
			float p[3];
			p[0] = (float)(vertex[i][0] + v[0]);
			p[1] = (float)(vertex[i][1] + v[1]);
			p[2] = (float)(vertex[i][2] + v[2]);
			mesh->barycentricCoor(ravine_point[ravine_N].id, p, ravine_point[ravine_N].u, ravine_point[ravine_N].v);
			ravine_N++;
		}
	}
}

void FeatureDetector::generateRidgeLine(float length, float step, double T1, double T2)
{
	int ridge_point_N = mesh->ridge_point_N;
	mesh->ridge_line = new MeshData::tri_point*[ridge_point_N];
	mesh->ridge_line_N = new int[ridge_point_N];
	int i;  for(i = 0; i<ridge_point_N; i++){
		int *f = mesh->face[mesh->ridge_point[i].id];
		double t1 = mesh->ridge_point[i].u;
		double t2 = mesh->ridge_point[i].v;
		double t0 = 1.0 - t1 - t2;
		double k = fabs(t0*mesh->k_max[f[0]] + t1*mesh->k_max[f[1]] + t2*mesh->k_max[f[2]]);
		mesh->traceRidgeDir(mesh->ridge_point[i], length*(float)k, mesh->ridge_line[i], mesh->ridge_line_N[i], step, T1);
	}

	int ravine_point_N = mesh->ravine_point_N;
	mesh->ravine_line = new MeshData::tri_point*[ravine_point_N];
	mesh->ravine_line_N = new int[ravine_point_N];
	for(i=0; i<ravine_point_N; i++){
		int *f = mesh->face[mesh->ridge_point[i].id];
		double t1 = mesh->ravine_point[i].u;
		double t2 = mesh->ravine_point[i].v;
		double t0 = 1.0 - t1 - t2;
		double k = fabs(t0*mesh->k_min[f[0]] + t1*mesh->k_min[f[1]] + t2*mesh->k_min[f[2]]);
		mesh->traceRavineDir(mesh->ravine_point[i], length*(float)k, mesh->ravine_line[i], mesh->ravine_line_N[i], step, T2);
	}
		
}

void FeatureDetector::generateRidgeLine2(float length, float step, double T1, double T2)
{
	int ridge_point_N = mesh->ridge_point_N;
	mesh->ridge_line = new MeshData::tri_point*[ridge_point_N];
	mesh->ridge_line_N = new int[ridge_point_N];
	int i;  for(i = 0; i<ridge_point_N; i++)
		mesh->traceRidgeDir2(mesh->ridge_point[i], length, mesh->ridge_line[i], mesh->ridge_line_N[i], step, T1);

	int ravine_point_N = mesh->ravine_point_N;
	mesh->ravine_line = new MeshData::tri_point*[ravine_point_N];
	mesh->ravine_line_N = new int[ravine_point_N];
	for(i=0; i<ravine_point_N; i++)
		mesh->traceRavineDir2(mesh->ravine_point[i], length, mesh->ravine_line[i], mesh->ravine_line_N[i], step, T2);
	
}

void FeatureDetector::automaticThresholding(float hi, float low, double &ridge_low, double &ravine_low)
{
	int vertex_N = mesh->vertex_N;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int *isBound = mesh->isBound;

	float *k1 = new float[vertex_N];
	float *k2 = new float[vertex_N];

	int c_max = 0;
	int c_min = 0;
	int i;  for(i = 0; i<vertex_N; i++){
		if(k_max[i] > 0)
			k1[c_max++] = (float)k_max[i];
		if(k_min[i] < 0)
			k2[c_min++] = (float)k_min[i];
	}
	quickSort(k1, 0, c_max-1);
	quickSort(k2, 0, c_min-1);

	float ridge_hi = k1[(int)(hi*(c_max-1))];
	ridge_low = k1[(int)(low*(c_max-1))];

	float ravine_hi = k2[(int)((1.0-hi)*(c_min-1))];
	ravine_low = k2[(int)((1.0-low)*(c_min-1))];

	delete[] k1;
	delete[] k2;

	for(i=0; i<vertex_N; i++){
		if(isBound[i])
			continue;

		if(isRidge[i]){
			if(k_max[i] > ridge_hi)
				isRidge[i] = true;
			else
				isRidge[i] = false;
		}

		if(isRavine[i]){
			if(k_min[i] < ravine_hi)
				isRavine[i] = true;
			else
				isRavine[i] = false;
		}
	}
}

void FeatureDetector::crossThresholding()
{
	int vertex_N = mesh->vertex_N;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int *isBound = mesh->isBound;

	int i;  for(i = 0; i<vertex_N; i++){
		if(isRidge[i]){
			if(fabs(k_max[i]) < fabs(k_min[i]))
				isRidge[i] = false;
		}
		if(isRavine[i]){
			if(fabs(k_max[i]) > fabs(k_min[i]))
				isRavine[i] = false;
		}
	}
}

void FeatureDetector::connectRR(float T, double T1, double T2, BOOL isCross)
{
	int vertex_N = mesh->vertex_N;
	int *degree = new int[vertex_N];
	int *id = new int[vertex_N];
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;

	int i;  for(i = 0; i<vertex_N; i++){
		if(isRidge[i]){
			degree[i] = 0;
			id[i] = 0;
			int *nei, nei_N;
			mesh->getConsistent1Ring(i, nei, nei_N);
			if(nei_N == 0)
				continue;
			int j;  for(i = 0; j<nei_N; j++)
				if(isRidge[nei[j]])
					degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}

	int current_id = 1;
	for(i=0; i<vertex_N; i++)
		if(id[i] == 0){
			labelConnectedComponent(i, 0, current_id, id);
			current_id++;
		}
			

	BOOL update = true;
	while(update){
		update = false;
		for(i=0; i<vertex_N; i++){
			if(degree[i] != 0 && degree[i] != 1)
				continue;
			Node* path = searchPathRidge(i, T, mesh->k_max, 5, degree, id, T1, isCross, mesh->k_min);
			if(path != NULL){
				int new_id = id[path->v];
				degree[path->v]++;
				for(Node* current = path->next; current->v!=i; current=current->next){
					isRidge[current->v] = true;
					degree[current->v] = 2;
					id[current->v] = new_id;
				}
				degree[i]++;
				labelConnectedComponent(i, id[i], new_id, id);
				update = true;
				delete path;
			}
		}
		update = false;
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRidge[i] = false;

	//////////////////Ravine///////////////////
	for(i=0; i<vertex_N; i++){
		if(isRavine[i]){
			degree[i] = 0;
			id[i] = 0;
			int *nei, nei_N;
			mesh->getConsistent1Ring(i, nei, nei_N);
			if(nei_N == 0)
				continue;
			int j;  for(i = 0; j<nei_N; j++)
				if(isRavine[nei[j]])
					degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -1;
		}
	}

	current_id = 1;
	for(i=0; i<vertex_N; i++)
		if(id[i] == 0){
			labelConnectedComponent(i, 0, current_id, id);
			current_id++;
		}
			

	update = true;
	while(update){
		update = false;
		for(i=0; i<vertex_N; i++){
			if(degree[i] != 0 && degree[i] != 1)
				continue;
			Node* path = searchPathRavine(i, T, mesh->k_min, -5, degree, id, T2, isCross, mesh->k_max);
			if(path != NULL){
				int new_id = id[path->v];
				degree[path->v]++;
				for(Node* current = path->next; current->v!=i; current=current->next){
					isRavine[current->v] = true;
					degree[current->v] = 2;
					id[current->v] = new_id;
				}
				degree[i]++;
				labelConnectedComponent(i, id[i], new_id, id);
				update = true;
				delete path;
			}
		}
		update = false;
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRavine[i] = false;

	delete[] id;
	delete[] degree;
}

Node* FeatureDetector::searchPathRidge(int start, float T, double *k, float K, int *degree, int *id, double KT, BOOL isCross, double* K2)
{
	int vertex_N = mesh->vertex_N;
	float (*vertex)[3] = mesh->vertex;
	float* label = new float[vertex_N];
	int* pred = new int[vertex_N];

	int* heap = new int[vertex_N];
	int* index = new int[vertex_N];
	int last_heap = 0;
	int i;  for(i = 0; i<vertex_N; i++)
		index[i] = -1;

	for(i=0; i<vertex_N; i++){
		label[i] = 1000000;
		pred[i] = i;
	}
	label[start] = 0;
	
	int u = start;
	int last_u = -1;
	int end = -1;
	while(u != last_u){
		int *nei, nei_N;
		mesh->getConsistent1Ring(u, nei, nei_N);
		for(i=0; i<nei_N; i++){
			int v = nei[i];
			if(k[v] < KT)
				continue;
			if(isCross && fabs(k[v]) < fabs(K2[v]))
				continue;
			if(id[start] == id[v])
				continue;

			float w_k = (float)(0.5*((k[u]-K)*(k[u]-K) + (k[v]-K)*(k[v]-K)) 
				+ (k[u]-k[v])*(k[u]-k[v])/6.0);// + 1.0;
			float w = (float)(w_k*MeshData::DIST(vertex[u], vertex[v]));
			
			float M = label[u] + w;
			if(index[v] < 0){
				label[v] = M;
				pred[v] = u;
				//insert;
				heap[++last_heap] = v;
				index[v] = last_heap;
				//upheap;
				upheap(label, last_heap, last_heap, heap, index);
			}
			else if(M < label[v]){
				label[v] = M;
				pred[v] = u;
				//change;
				if(index[v] != 1 && M < label[heap[index[v]/2]])
					//upheap;
					upheap(label, last_heap, index[v], heap, index);
				else
					//downheap;
					downheap(label, last_heap, index[v], heap, index);
			}
		}

		if(u == end || last_heap == 0 )
			break;		

		//delete;
		last_u = u;
		u = heap[1];
		if((degree[u] == 0 || degree[u] == 1)  && id[u] != id[start]){
			end = u;
			break;
		}
		else if(label[u] > T)
			break;
		else if(last_heap == 0)
			break;
		heap[1] = heap[last_heap--];
		index[heap[1]] = 1;
		//down heap;
		downheap(label, last_heap, 1, heap, index);
	}
	delete[] index;
	delete[] heap;
	delete[] label;

	
	Node* path = NULL;
	if(end != -1){
		path = new Node;
		for(int current = end; current != start; current = pred[current])
			path->append(current, -1);
		path->append(start, -1);
	}

	delete[] pred;

	return path;
}

Node* FeatureDetector::searchPathRavine(int start, float T, double *k, float K, int *degree, int *id, double KT, BOOL isCross, double *K2)
{
	int vertex_N = mesh->vertex_N;
	float (*vertex)[3] = mesh->vertex;
	float* label = new float[vertex_N];
	int* pred = new int[vertex_N];

	int* heap = new int[vertex_N];
	int* index = new int[vertex_N];
	int last_heap = 0;
	int i;  for(i = 0; i<vertex_N; i++)
		index[i] = -1;

	for(i=0; i<vertex_N; i++){
		label[i] = 1000000;
		pred[i] = i;
	}
	label[start] = 0;
	
	int u = start;
	int last_u = -1;
	int end = -1;
	while(u != last_u){
		int *nei, nei_N;
		mesh->getConsistent1Ring(u, nei, nei_N);
		for(i=0; i<nei_N; i++){
			int v = nei[i];
			if(k[v] > KT)
				continue;
			if(isCross && fabs(k[v]) < fabs(K2[v]))
				continue;
			if(id[start] == id[v])
				continue;

			float w_k = (float)(0.5*((k[u]-K)*(k[u]-K) + (k[v]-K)*(k[v]-K)) 
				+ (k[u]-k[v])*(k[u]-k[v])/6.0);// + 1.0;
			float w = (float)(w_k*MeshData::DIST(vertex[u], vertex[v]));
			
			float M = label[u] + w;
			if(index[v] < 0){
				label[v] = M;
				pred[v] = u;
				//insert;
				heap[++last_heap] = v;
				index[v] = last_heap;
				//upheap;
				upheap(label, last_heap, last_heap, heap, index);
			}
			else if(M < label[v]){
				label[v] = M;
				pred[v] = u;
				//change;
				if(index[v] != 1 && M < label[heap[index[v]/2]])
					//upheap;
					upheap(label, last_heap, index[v], heap, index);
				else
					//downheap;
					downheap(label, last_heap, index[v], heap, index);
			}
		}

		if(u == end || last_heap == 0)
			break;		

		//delete;
		last_u = u;
		u = heap[1];
		if((degree[u] == 0 || degree[u] == 1)  && id[u] != id[start]){
			end = u;
			break;
		}
		else if(label[u] > T)
			break;
		else if(last_heap == 0)
			break;
		heap[1] = heap[last_heap--];
		index[heap[1]] = 1;
		//down heap;
		downheap(label, last_heap, 1, heap, index);
	}
	delete[] index;
	delete[] heap;
	delete[] label;

	
	Node* path = NULL;
	if(end != -1){
		path = new Node;
		for(int current = end; current != start; current = pred[current])
			path->append(current, -1);
		path->append(start, -1);
	}

	delete[] pred;

	return path;
}

void FeatureDetector::addBridgeRR(double T1, double T2)
{
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int vertex_N = mesh->vertex_N;
	int **link = mesh->vertex_link_v;
	int *degree = mesh->degree_v;

	//Ridge
	for(int i = 0; i<vertex_N; i++){
		if(isRidge[i])
			continue;

		int* l = link[i];
		int deg = degree[i];
		for(int i = 0; j<deg; j++){
			if(isRidge[l[j]]){
				int pair = l[(j+1)%deg];
				int k = l[(j+2)%deg];
				if(!isRidge[pair] && isRidge[k]){
					if(fabs(mesh->k_max[i]) > fabs(mesh->k_max[pair]))
						isRidge[i] = true;
					else
						isRidge[pair] = true;
				}
			}
		}
	}

	//Ravine
	for(int i=0; i<vertex_N; i++){
		if(isRavine[i])
			continue;

		int* l = link[i];
		int deg = degree[i];
		for(int i = 0; j<deg; j++){
			if(isRavine[l[j]]){
				int pair = l[(j+1)%deg];
				int k = l[(j+2)%deg];
				if(!isRavine[pair] && isRavine[k]){
					if(fabs(mesh->k_min[i]) > fabs(mesh->k_min[pair]))
						isRavine[i] = true;
					else
						isRavine[pair] = true;
				}
			}
		}
	}
}

void FeatureDetector::setRidgeEdges()
{
	int vertex_N = mesh->vertex_N;
	int (*face)[3] = mesh->face;
	int face_N = mesh->face_N;

	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;

	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;

	Node** ridge_edge = mesh->ridge_edge;
	Node** ravine_edge = mesh->ravine_edge;

	if(ridge_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ridge_edge[i];
		delete[] ridge_edge;
	}
	
	if(ravine_edge != NULL){
		int i;  for(i = 0; i<vertex_N; i++)
			delete ravine_edge[i];
		delete[] ravine_edge;
	}

	ridge_edge = mesh->ridge_edge = new Node*[vertex_N];
	ravine_edge = mesh->ravine_edge = new Node*[vertex_N];

	int i;  for(i = 0; i<vertex_N; i++){
		ridge_edge[i] = new Node;
		ravine_edge[i] = new Node;
	}

	int ridge_N = 0;
	int ravine_N = 0;
	for(i=0; i<vertex_N; i++){
		if(isRidge[i])
			ridge_N++;
		if(isRavine[i])
			ravine_N++;
	}


	//Ridge
	int *index = new int[ridge_N];
	float *tmp_k = new float[ridge_N];
	ridge_N = 0;
	for(i=0; i<vertex_N; i++){
		if(isRidge[i]){
			index[ridge_N] = i;
			tmp_k[ridge_N] = -(float)k_max[i];
			ridge_N++;
		}
	}
	quickSort(index, tmp_k, 0, ridge_N-1);
	delete[] tmp_k;

	for(i=0; i<ridge_N; i++){
		int j = index[i];
		if(ridge_edge[j]->next != NULL)
			continue;
		double K = k_max[j];
		connectRidgeBF(j, K);
	}

	BOOL* isVisit = new BOOL[vertex_N];
	
	for(i=0; i<vertex_N; i++)
		isVisit[i] = false;
	for(i=0; i<ridge_N; i++){
		int j = index[i];
		if(isVisit[j])
			continue;
		double K = k_max[j];
		checkBranchRidge(j, K, isVisit, 2.5);
	}

	delete[] index;

	//Ravine
	index = new int[ravine_N];
	tmp_k = new float[ravine_N];
	ravine_N = 0;
	for(i=0; i<vertex_N; i++){
		if(isRavine[i]){
			index[ravine_N] = i;
			tmp_k[ravine_N] = (float)k_min[i];
			ravine_N++;
		}
	}
	quickSort(index, tmp_k, 0, ravine_N-1);
	delete[] tmp_k;

	for(i=0; i<ravine_N; i++){
		int j = index[i];
		if(ravine_edge[j]->next != NULL)
			continue;
		double K = k_min[j];
		connectRavineBF(j, K);
	}

	for(i=0; i<vertex_N; i++)
		isVisit[i] = false;
	for(i=0; i<ravine_N; i++){
		int j = index[i];
		if(isVisit[j])
			continue;
		double K = k_min[j];
		checkBranchRavine(j, K, isVisit, 2.5);
	}

	delete[] index;
	delete[] isVisit;
/*
	for(i=0; i<face_N; i++){
		int *f = face[i];
		int skip = -1;

		if(isRidge[f[0]] && isRidge[f[1]] && isRidge[f[2]]){
			double e1 = k_max[f[0]] + k_max[f[1]];
			double e2 = k_max[f[1]] + k_max[f[2]];
			double e3 = k_max[f[2]] + k_max[f[0]];

			if(e1 < e2 && e1 < e3)
				skip = 0;
			else if(e2 < e3 && e2 < e1)
				skip = 1;
			else
				skip = 2;
		}
		int j;  for(i = 0; j<3; j++){
			if(j == skip)
				continue;
			int i1 = f[j];
			int i2 = f[(j+1)%3];
			if(isRidge[i1] && isRidge[i2]){
				ridge_edge[i1]->append(i2, -1);
				//ridge_edge[i2]->append(i1, -1);
			}
		}
	}

	//Ravine
	for(i=0; i<face_N; i++){
		int *f = face[i];
		int skip = -1;

		if(isRavine[f[0]] && isRavine[f[1]] && isRavine[f[2]]){
			double e1 = k_min[f[0]] + k_min[f[1]];
			double e2 = k_min[f[1]] + k_min[f[2]];
			double e3 = k_min[f[2]] + k_min[f[0]];

			if(e1 > e2 && e1 > e3)
				skip = 0;
			else if(e2 > e3 && e2 > e1)
				skip = 1;
			else
				skip = 2;
		}
		int j;  for(i = 0; j<3; j++){
			if(j == skip)
				continue;
			int i1 = f[j];
			int i2 = f[(j+1)%3];
			if(isRavine[i1] && isRavine[i2]){
				ravine_edge[i1]->append(i2, -1);
				//ravine_edge[i2]->append(i1, -1);
			}
		}
	}*/
}

void FeatureDetector::quickSortD(int *index, double *w, int start, int end)
{
	if(start < end){
		double v = w[end];
		int i = start-1;
		int j = end;
		while(j > i){
			for(i = i+1; w[i] < v; i++);
			for(j = j-1; w[j] > v; j--);
			double t = w[i]; 
			w[i] = w[j]; 
			w[j] = t;

			int tmp = index[i];
			index[i] = index[j];
			index[j] = tmp;
		}
		double t = w[j];
		w[j] = w[i];
		w[i] = w[end];
		w[end] = t;

		int tmp = index[j];
		index[j] = index[i];
		index[i] = index[end];
		index[end] = tmp;

		quickSortD(index, w, start, i-1);
		quickSortD(index, w, i+1, end);
	}
	else
		return;
}

double FeatureDetector::connectRidgeDF(int i, double K)
{
	BOOL* isRidge = mesh->isRidge;
	int* link = mesh->vertex_link_v[i];
	int deg = mesh->degree_v[i];
	Node** ridge_edge = mesh->ridge_edge;

	int j;  for(i = 0; j<deg; j++){
		int k = link[j];
		if(isRidge[k] && ridge_edge[k]->next == NULL){
			ridge_edge[i]->append(k, -1);
			ridge_edge[k]->append(i, -1);
			connectRidgeDF(k, K);
		}
	}
	return 0;
}

double FeatureDetector::connectRavineDF(int i, double K)
{
	BOOL* isRavine = mesh->isRavine;
	int* link = mesh->vertex_link_v[i];
	int deg = mesh->degree_v[i];
	Node** ravine_edge = mesh->ravine_edge;

	int j;  for(i = 0; j<deg; j++){
		int k = link[j];
		if(isRavine[k] && ravine_edge[k]->next == NULL){
			ravine_edge[i]->append(k, -1);
			ravine_edge[k]->append(i, -1);
			connectRavineDF(k, K);
		}
	}
	return 0;
}

void FeatureDetector::connectRidgeBF(int i, double K)
{
	BOOL* isRidge = mesh->isRidge;
	int** link = mesh->vertex_link_v;
	int* degree = mesh->degree_v;
	Node** ridge_edge = mesh->ridge_edge;
	Node* Q = new Node;
	Q->append(i, -1);
	double* k_max = mesh->k_max;
	float (*vertex)[3] = mesh->vertex;

	while(Q->next != NULL){
		int index = Q->v;
		Node* tmp = Q;
		Q->next->tail = Q->tail;
		Q = Q->next;
		tmp->next = NULL;
		delete tmp;

		int deg = degree[index];
		int* l = link[index];
		int count = 0;
		double* w = new double[deg];
		int* id = new int[deg];
		int j;  for(i = 0; j<deg; j++){
			int k = l[j];
			if(isRidge[k] && ridge_edge[k]->next == NULL){
				w[count] = 0.5*((k_max[index]-K)*(k_max[index]-K) 
									+ (k_max[k]-K)*(k_max[k]-K)) 
							+ (k_max[index]-k_max[k])*(k_max[index]-k_max[k])/6.0;
				w[count] *= MeshData::DIST(vertex[k], vertex[index]);
				id[count] = k;
				count++;

				ridge_edge[index]->append(k, -1);
				ridge_edge[k]->append(index, -1);
			}
		}
		for(int m=0; m<count; m++){
			int min = 0;
			for(int n=1; n<count; n++){
				if(w[min] > w[n])
					min = n;
			}
			w[min] = 100000000;
			Q->append(id[min], -1);
		}

		delete[] w;
		delete[] id;
	}

	delete Q;
}

void FeatureDetector::connectRavineBF(int i, double K)
{
	BOOL* isRavine = mesh->isRavine;
	int** link = mesh->vertex_link_v;
	int* degree = mesh->degree_v;
	Node** ravine_edge = mesh->ravine_edge;
	Node* Q = new Node;
	Q->append(i, -1);
	double* k_min = mesh->k_min;
	float (*vertex)[3] = mesh->vertex;

	while(Q->next != NULL){
		int index = Q->v;
		Node* tmp = Q;
		Q->next->tail = Q->tail;
		Q = Q->next;
		tmp->next = NULL;
		delete tmp;

		int deg = degree[index];
		int* l = link[index];
		int count = 0;
		double* w = new double[deg];
		int* id = new int[deg];
		int j;  for(i = 0; j<deg; j++){
			int k = l[j];
			if(isRavine[k] && ravine_edge[k]->next == NULL){
				w[count] = 0.5*((k_min[index]-K)*(k_min[index]-K) 
									+ (k_min[k]-K)*(k_min[k]-K)) 
							+ (k_min[index]-k_min[k])*(k_min[index]-k_min[k])/6.0;
				w[count] *= MeshData::DIST(vertex[k], vertex[index]);
				id[count] = k;
				count++;

				ravine_edge[index]->append(k, -1);
				ravine_edge[k]->append(index, -1);
			}
		}
		for(int m=0; m<count; m++){
			int min = 0;
			for(int n=1; n<count; n++){
				if(w[min] > w[n])
					min = n;
			}
			w[min] = 100000000;
			Q->append(id[min], -1);
		}

		delete[] w;
		delete[] id;
	}

	delete Q;
}

double FeatureDetector::checkBranchRidge(int index, double K, BOOL *isVisit, double T)
{
	Node** ridge_edge = mesh->ridge_edge;
	float (*vertex)[3] = mesh->vertex;
	double* k_max = mesh->k_max;
	int deg = mesh->degree_v[index];
	isVisit[index] = true;

	double total_w = 0;
	int count = 0;
	int* link = new int[deg];
	double* w = new double[deg];
	for(Node* current = ridge_edge[index]; current->next != NULL; current = current->next){
		int i = current->v;
		link[count] = i;
		w[count] = 100000000;
		if(isVisit[i]){
			count++;
			continue;
		}
		/*
		w[count] = 0.5*((k_max[index]-K)*(k_max[index]-K) 
									+ (k_max[i]-K)*(k_max[i]-K)) 
					+ (k_max[index]-k_max[i])*(k_max[index]-k_max[i])/6.0;*/
		w[count] = 1; //MeshData::DIST(vertex[i], vertex[index]);
		w[count] += checkBranchRidge(i, K, isVisit, T);
		total_w += w[count];
		
		count++;
	}
	if(count > 2){
		int j;  for(i = 0; j<count-2; j++){
		int min = 1;
		for(int i=1; i<count; i++){
			if(w[min] > w[i])
				min = i;
		}
		if(w[min] < T)
			eliminateBranchRidge(index, link[min]);
		total_w -= w[min];
		w[min] = 1000000;
		}
	}

	delete[] link;
	delete[] w;

	return total_w;
}

void FeatureDetector::eliminateBranchRidge(int from, int to)
{
	Node** ridge_edge = mesh->ridge_edge;

	Node* tmp = new Node;
	for(Node* current=ridge_edge[to]; current->next != NULL; current=current->next)
		tmp->append(current->v, -1);

	for(current=tmp; current->next != NULL; current=current->next)
		if(current->v != from)
			eliminateBranchRidge(to, current->v);

	delete tmp;

	tmp = new Node;
	for(current=ridge_edge[to]; current->next != NULL; current=current->next){
		if(current->v != from)
			tmp->append(current->v, -1);
	}
	delete ridge_edge[to];
	ridge_edge[to] = tmp;

	tmp = new Node;
	for(current=ridge_edge[from]; current->next != NULL; current=current->next){
		if(current->v != to)
			tmp->append(current->v, -1);
	}
	delete ridge_edge[from];
	ridge_edge[from] = tmp;
}

double FeatureDetector::checkBranchRavine(int index, double K, BOOL *isVisit, double T)
{
	Node** ravine_edge = mesh->ravine_edge;
	float (*vertex)[3] = mesh->vertex;
	double* k_min = mesh->k_min;
	int deg = mesh->degree_v[index];
	isVisit[index] = true;

	double total_w = 0;
	int count = 0;
	int* link = new int[deg];
	double* w = new double[deg];
	for(Node* current = ravine_edge[index]; current->next != NULL; current = current->next){
		int i = current->v;
		link[count] = i;
		w[count] = 100000000;
		if(isVisit[i]){
			count++;
			continue;
		}
		/*
		w[count] = 0.5*((k_min[index]-K)*(k_min[index]-K) 
									+ (k_min[i]-K)*(k_min[i]-K)) 
					+ (k_min[index]-k_min[i])*(k_min[index]-k_min[i])/6.0;*/
		w[count] = 1; //MeshData::DIST(vertex[i], vertex[index]);
		w[count] += checkBranchRavine(i, K, isVisit, T);
		total_w += w[count];
		
		count++;
	}
	if(count > 2){
		int j;  for(i = 0; j<count-2; j++){
		int min = 1;
		for(int i=1; i<count; i++){
			if(w[min] > w[i])
				min = i;
		}
		if(w[min] < T)
			eliminateBranchRavine(index, link[min]);
		total_w -= w[min];
		w[min] = 1000000;
		}
	}

	delete[] link;
	delete[] w;

	return total_w;
}

void FeatureDetector::eliminateBranchRavine(int from, int to)
{
	Node** ravine_edge = mesh->ravine_edge;

	Node* tmp = new Node;
	for(Node* current=ravine_edge[to]; current->next != NULL; current=current->next)
		tmp->append(current->v, -1);

	for(current=tmp; current->next != NULL; current=current->next)
		if(current->v != from)
			eliminateBranchRavine(to, current->v);

	delete tmp;

	tmp = new Node;
	for(current=ravine_edge[to]; current->next != NULL; current=current->next){
		if(current->v != from)
			tmp->append(current->v, -1);
	}
	delete ravine_edge[to];
	ravine_edge[to] = tmp;

	tmp = new Node;
	for(current=ravine_edge[from]; current->next != NULL; current=current->next){
		if(current->v != to)
			tmp->append(current->v, -1);
	}
	delete ravine_edge[from];
	ravine_edge[from] = tmp;
}

void FeatureDetector::connectRRedge(float T, double T_ridge, double T_ravine)
{
	int vertex_N = mesh->vertex_N;
	int *degree = new int[vertex_N];
	int *id = new int[vertex_N];
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	Node **ridge_edge = mesh->ridge_edge;
	Node **ravine_edge = mesh->ravine_edge;
	float *strength;
	int *index = new int[vertex_N];
	float *end_strength = new float[vertex_N];
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	float (*vertex)[3] = mesh->vertex;

	int i;  for(i = 0; i<vertex_N; i++){
		if(ridge_edge[i]->next == NULL)
			isRidge[i] = false;
		if(ravine_edge[i]->next == NULL)
			isRavine[i] = false;
	}

	for(i=0; i<vertex_N; i++){
		if(isRidge[i]){
			degree[i] = 0;
			id[i] = -1;
			for(Node* current=ridge_edge[i]; current->next!=NULL; current=current->next)
				degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -2;
		}
	}

	int current_id = 0;
	for(i=0; i<vertex_N; i++)
		if(id[i] == -1){
			labelConnectedComponent(i, -1, current_id, id);
			current_id++;
		}

	strength = new float[current_id];
	for(i=0; i<current_id; i++)
		strength[i] = 0;
	for(i=0; i<vertex_N; i++){
		if(id[i] < 0)
			continue;
		if(degree[i] == 0){
			strength[id[i]] = (float)(1.0/fabs(k_max[i]));
		}
		else{
			for(Node* current=ridge_edge[i]; current->next!=NULL; current=current->next)
				strength[id[i]] -= (float)(0.5*fabs(k_max[i])*MeshData::DIST(vertex[i], vertex[current->v]));
		}
	}
	int end_N = 0;
	for(i=0; i<vertex_N; i++){
		if(degree[i] == 0 || degree[i] == 1){
			end_strength[end_N] = strength[id[i]];
			index[end_N] = i;
			end_N++;
		}
	}
	quickSort(index, end_strength, 0, end_N-1);

	int j;  for(i = 0; j<end_N; j++){
		i = index[j];
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPathRidge(i, T, mesh->k_max, 30, degree, id, T_ridge, true, mesh->k_min);
		if(path != NULL){
			int new_id = id[path->v];
			degree[path->v]++;
			ridge_edge[path->v]->append(path->next->v, -1);
			ridge_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRidge[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ridge_edge[current->v]->append(current->next->v, -1);
				ridge_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
	}
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRidge[i] = false;

	delete[] strength;

	//////////////////Ravine///////////////////
	for(i=0; i<vertex_N; i++){
		if(isRavine[i]){
			degree[i] = 0;
			id[i] = -1;
			for(Node* current=ravine_edge[i]; current->next!=NULL; current=current->next)
				degree[i]++;
		}
		else{
			degree[i] = -1;
			id[i] = -2;
		}
	}

	current_id = 0;
	for(i=0; i<vertex_N; i++)
		if(id[i] == -1){
			labelConnectedComponent(i, -1, current_id, id);
			current_id++;
		}

	strength = new float[current_id];
	for(i=0; i<current_id; i++)
		strength[i] = 0;
	for(i=0; i<vertex_N; i++){
		if(id[i] < 0)
			continue;
		if(degree[i] == 0){
			strength[id[i]] = (float)(1.0/fabs(mesh->k_min[i]));
		}
		else{
			for(Node* current=ravine_edge[i]; current->next!=NULL; current=current->next)
				strength[id[i]] -= (float)(0.5*fabs(mesh->k_min[i])*MeshData::DIST(vertex[i], vertex[current->v]));
		}
	}
	end_N = 0;
	for(i=0; i<vertex_N; i++){
		if(degree[i] == 0 || degree[i] == 1){
			end_strength[end_N] = strength[id[i]];
			index[end_N] = i;
			end_N++;
		}
	}
	quickSort(index, end_strength, 0, end_N-1);
			
	for(j=0; j<end_N; j++){
		i = index[j];
		if(degree[i] != 0 && degree[i] != 1)
			continue;
		Node* path = searchPathRavine(i, T, mesh->k_min, -30, degree, id, T_ravine, true, mesh->k_max);
		if(path != NULL){
			int new_id = id[path->v];
			degree[path->v]++;
			ravine_edge[path->v]->append(path->next->v, -1);
			ravine_edge[path->next->v]->append(path->v, -1);
			for(Node* current = path->next; current->v!=i; current=current->next){
				isRavine[current->v] = true;
				degree[current->v] = 2;
				id[current->v] = new_id;
				ravine_edge[current->v]->append(current->next->v, -1);
				ravine_edge[current->next->v]->append(current->v, -1);
			}
			degree[i]++;
			labelConnectedComponent(i, id[i], new_id, id);
			delete path;
		}
	}
	
	for(i=0; i<vertex_N; i++)
		if(degree[i] == 0)
			isRavine[i] = false;

	delete[] id;
	delete[] degree;
	delete[] index;
	delete[] strength;
	delete[] end_strength;
}

void FeatureDetector::ridgeTriangle()
{
	int vertex_N = mesh->vertex_N;
	int face_N = mesh->face_N;
	BOOL *ridge_tri = mesh->ridge_tri = new BOOL[face_N];
	BOOL *ravine_tri = mesh->ravine_tri = new BOOL[face_N];
	BOOL *isRidge = mesh->isRidge;
	BOOL *isRavine = mesh->isRavine;
	int *degree = mesh->degree_f;
	int **link = mesh->vertex_link_f;
	int **link_v = mesh->vertex_link_v;

	int i;  for(i = 0; i<face_N; i++){
		ridge_tri[i] = false;
		ravine_tri[i] = false;
	}

	for(i=0; i<vertex_N; i++){
		int *l = link[i];
		if(isRidge[i]){
			int j;  for(i = 0; j<degree[i]; j++){
				if(isRidge[link_v[i][j]] || isRidge[link_v[i][(j+1)%degree[i]]])
					ridge_tri[l[j]] = true;
			}
		}
		if(isRavine[i]){
			int j;  for(i = 0; j<degree[i]; j++)
				if(isRavine[link_v[i][j]] || isRavine[link_v[i][(j+1)%degree[i]]])
					ravine_tri[l[j]] = true;
		}
	}

}

void FeatureDetector::ridgeBelyaevCGI98()
{
	double* Kmax = mesh->k_max;
	double* Kmin = mesh->k_min;
	double* Rmax = mesh->r_max;
	double* Rmin = mesh->r_min;
	double* Dmax = mesh->d_max;
	double* Dmin = mesh->d_min;
	double (*Tmax)[3] = mesh->t_max;
	double (*Tmin)[3] = mesh->t_min;

	float (*vertex)[3] = mesh->vertex;
	int (*face)[3] = mesh->face;
	int faceN = mesh->face_N;

	RList* ridge = new RList;
	RList* ravine = new RList;
	ridge->next = NULL;
	ravine->next = NULL;

	int i;  for(i = 0; i<faceN; i++){
		int* f = face[i];
	
		if(MeshData::LENGTH(Tmax[f[0]]) == 0 || MeshData::LENGTH(Tmax[f[1]]) == 0 || MeshData::LENGTH(Tmax[f[2]]) == 0)
			continue;

		bool first1 = false;
		bool first2 = false;
		bool second1 = false;
		bool second2 = false;
		RList* list1;
		RList* list2;
		int j;  for(i = 0; j<3; j++){
			int i1 = f[j];
			int i2 = f[(j+1)%3];

			//double tmp = Rmax[i2];
			if(MeshData::DOT(Tmax[i1], Tmax[i2]) < 0){
				Rmax[i2] = -Rmax[i2];
				Tmax[i2][0] = -Tmax[i2][0];
				Tmax[i2][1] = -Tmax[i2][1];
				Tmax[i2][2] = -Tmax[i2][2];
			}
			if(Rmax[i1]*Rmax[i2] < 0 && Kmax[i1] > fabs(Kmin[i1]) && Kmax[i2] > fabs(Kmin[i2])){
				double w1 = fabs(Rmax[i2]);
				double w2 = fabs(Rmax[i1]);
				double w = w1 + w2;
				w1 /= w;
				w2 /= w;

				float x = (float)(w1*vertex[i1][0] + w2*vertex[i2][0]);
				float y = (float)(w1*vertex[i1][1] + w2*vertex[i2][1]);
				float z = (float)(w1*vertex[i1][2] + w2*vertex[i2][2]);
				double k = w1*Kmax[i1] + w2*Kmax[i2];
				double d; // = w1*Dmax[i1] + w2*Dmax[i2];
				
				float v[3];
				MeshData::VEC(v, vertex[i1], vertex[i2]);
				if(Rmax[i1]*MeshData::DOT(v, Tmax[i1]) 
					- Rmax[i2]*MeshData::DOT(v, Tmax[i2]) > 0)
					d = -1;
				else
					d = 1;
				/*
				if(MeshData::DOT(v, Tmax[i1]) > 0){
					if(Rmax[i1] > 0)
						d = -1;
					else
						d = 1;
				}
				else{
					if(Rmax[i1] > 0)
						d = 1;
					else
						d = -1;
				}*/

				if(!first1){
					list1 = new RList;
					first1 = true;
					list1->p1[0] = x;
					list1->p1[1] = y;
					list1->p1[2] = z;
					list1->k1 = k;
					list1->d1 = d;
				}
				else{
					second1 = true;
					list1->p2[0] = x;
					list1->p2[1] = y;
					list1->p2[2] = z;
					list1->k2 = k;
					list1->d2 = d;

					list1->vis = (list1->k1 > 0) && (list1->k2 > 0) && ((list1->d1 < 0) && (list1->d2 < 0));
				}
			}

			//tmp = Rmin[i2];
			//if(MeshData::DOT(Tmin[i1], Tmin[i2]) < 0)
				//tmp = -tmp;

			if(MeshData::DOT(Tmin[i1], Tmin[i2]) < 0){
				Rmin[i2] = -Rmin[i2];
				Tmin[i2][0] = -Tmin[i2][0];
				Tmin[i2][1] = -Tmin[i2][1];
				Tmin[i2][2] = -Tmin[i2][2];
			}
			if(Rmin[i1]*Rmin[i2] < 0 && -Kmin[i1] > fabs(Kmax[i1]) && -Kmin[i2] > fabs(Kmax[i2])){
				double w1 = fabs(Rmin[i2]);
				double w2 = fabs(Rmin[i1]);
				double w = w1 + w2;
				w1 /= w;
				w2 /= w;

				float x = (float)(w1*vertex[i1][0] + w2*vertex[i2][0]);
				float y = (float)(w1*vertex[i1][1] + w2*vertex[i2][1]);
				float z = (float)(w1*vertex[i1][2] + w2*vertex[i2][2]);
				double k = w1*Kmin[i1] + w2*Kmin[i2];
				double d;// = w1*Dmin[i1] + w2*Dmin[i2];
				
				float v[3];
				MeshData::VEC(v, vertex[i1], vertex[i2]);
				if(Rmin[i1]*MeshData::DOT(v, Tmin[i1]) 
					- Rmin[i2]*MeshData::DOT(v, Tmin[i2]) < 0)
					d = 1;
				else
					d = -1;

				/*
				if(MeshData::DOT(v, Tmin[i1]) > 0){
					if(Rmin[i1] > 0)
						d = -1;
					else
						d = 1;
				}
				else{
					if(Rmin[i1] > 0)
						d = 1;
					else
						d = -1;
				}*/

				if(!first2){
					list2 = new RList;
					first2 = true;
					list2->p1[0] = x;
					list2->p1[1] = y;
					list2->p1[2] = z;
					list2->k1 = k;
					list2->d1 = d;
				}
				else{
					second2 = true;
					list2->p2[0] = x;
					list2->p2[1] = y;
					list2->p2[2] = z;
					list2->k2 = k;
					list2->d2 = d;

					list2->vis = (list2->k1 < 0) && (list2->k2 < 0) && ((list2->d1 > 0) && (list2->d2 > 0));
				}
			}
		}
		if(second1){
			list1->next = ridge;
			ridge = list1;
		}
		else if(first1){
			delete list1;
		}

		if(second2){
			list2->next = ravine;
			ravine = list2;
		}
		else if(first2){
			delete list2;
		}
	}

	mesh->ridge_L = ridge;
	mesh->ravine_L = ravine;
}

void FeatureDetector::ridgeBelyaevCGI98S()
{
	if(mesh->ridge_S != NULL)
		return;

	double* Kmax = mesh->k_max;
	double* Kmin = mesh->k_min;
	double* Rmax = mesh->r_max;
	double* Rmin = mesh->r_min;
	double* Dmax = mesh->d_max;
	double* Dmin = mesh->d_min;
	double (*Tmax)[3] = mesh->t_max;
	double (*Tmin)[3] = mesh->t_min;

	float (*vertex)[3] = mesh->vertex;
	int (*face)[3] = mesh->face;
	int faceN = mesh->face_N;

	int (*face_link)[3] = mesh->face_link_E;

	RPoint *(*ridgeP)[3] = new RPoint*[faceN][3];
	RPoint *(*ravineP)[3] = new RPoint*[faceN][3];

	int i;

	for(i=0; i<faceN; i++){
		int j;  for(i = 0; j<3; j++){
			ridgeP[i][j] = ravineP[i][j] = NULL;
		}
	}

	for(i=0; i<faceN; i++){
		int j;  for(i = 0; j<3; j++){
			int ad = face_link[i][j];
			if(i > ad)
				continue;
			int ad_i;
			if(face_link[ad][0] == i)
				ad_i = 0;
			else if(face_link[ad][1] == i)
				ad_i = 1;
			else
				ad_i = 2;

			int i1 = face[i][j];
			int i2 = face[i][(j+1)%3];
			float v[3];
			MeshData::VEC(v, vertex[i1], vertex[i2]);

			//ridge
			if(MeshData::DOT(Tmax[i1], Tmax[i2]) < 0){
				Rmax[i2] = -Rmax[i2];
				Tmax[i2][0] = -Tmax[i2][0];
				Tmax[i2][1] = -Tmax[i2][1];
				Tmax[i2][2] = -Tmax[i2][2];
			}

			if(//Kmax[i1]+Kmin[i1] > 0 && Kmax[i2]+Kmin[i2] > 0 &&
				Rmax[i1]*Rmax[i2] < 0 && 
				
				Kmax[i1] > fabs(Kmin[i1]) && 
				Kmax[i2] > fabs(Kmin[i2]) &&
				(Rmax[i1]*MeshData::DOT(v, Tmax[i1]) > 0 || 
				 - Rmax[i2]*MeshData::DOT(v, Tmax[i2]) > 0)){
				/*
	            Kmax[i1] > 0 && 
				Kmax[i2] > 0 &&
				(Dmax[i1] < 0 || Dmax[i2] < 0) 
				){*/

				RPoint* r = ridgeP[i][j] = ridgeP[ad][ad_i] = new RPoint;

				double w1 = fabs(Rmax[i2]);
				double w2 = fabs(Rmax[i1]);
				double w = w1 + w2;
				w1 /= w;
				w2 /= w;

				r->p[0] = (float)(w1*vertex[i1][0] + w2*vertex[i2][0]);
				r->p[1] = (float)(w1*vertex[i1][1] + w2*vertex[i2][1]);
				r->p[2] = (float)(w1*vertex[i1][2] + w2*vertex[i2][2]);
				r->k = w1*Kmax[i1] + w2*Kmax[i2];
				r->k2 = w1*Kmin[i1] + w2*Kmin[i2];
				double index = 2.0*atan((r->k + r->k2)/(r->k - r->k2))/PI;
				r->strong = true; //fabs(r->k) > fabs(r->k2); //index > 0; //(index > 0.25) && (index < 0.75);
				//r->strong = (Kmax[i1] > 2.0*fabs(Kmin[i1])) && (Kmax[i2] > 2.0*fabs(Kmin[i2]));
				//r->strong = r->k > fabs(k2);
				//r->k = fabs(r->k) - fabs(k2);
			}

			//ravine
			if(MeshData::DOT(Tmin[i1], Tmin[i2]) < 0){
				Rmin[i2] = -Rmin[i2];
				Tmin[i2][0] = -Tmin[i2][0];
				Tmin[i2][1] = -Tmin[i2][1];
				Tmin[i2][2] = -Tmin[i2][2];
			}

			if(//Kmax[i1]+Kmin[i1] < 0 && Kmax[i2]+Kmin[i2] < 0 &&
				Rmin[i1]*Rmin[i2] < 0 && 
				
				Kmin[i1] < -fabs(Kmax[i1]) && 
				Kmin[i2] < -fabs(Kmax[i2]) &&
				(Rmin[i1]*MeshData::DOT(v, Tmin[i1]) < 0 ||
				- Rmin[i2]*MeshData::DOT(v, Tmin[i2]) < 0)){
				/*
				Kmin[i1] < 0 && 
				Kmin[i2] < 0 && 
				(Dmin[i1] > 0 || Dmin[i2] > 0)
				){*/

				RPoint* r = ravineP[i][j] = ravineP[ad][ad_i] = new RPoint;

				double w1 = fabs(Rmin[i2]);
				double w2 = fabs(Rmin[i1]);
				double w = w1 + w2;
				w1 /= w;
				w2 /= w;

				r->p[0] = (float)(w1*vertex[i1][0] + w2*vertex[i2][0]);
				r->p[1] = (float)(w1*vertex[i1][1] + w2*vertex[i2][1]);
				r->p[2] = (float)(w1*vertex[i1][2] + w2*vertex[i2][2]);
				r->k = w1*Kmin[i1] + w2*Kmin[i2];
				r->k2 = w1*Kmax[i1] + w2*Kmax[i2];
				double index = 2.0*atan((r->k + r->k2)/(r->k - r->k2))/PI;
				r->strong =  true; //fabs(r->k) > fabs(r->k2); //index > 0; (index > 0.25) && (index < 0.75);
				//r->strong = 2.0*atan((r->k + r->k2)/(r->k - r->k2))/PI > 0.5;
				//r->strong = (-Kmin[i1] > 2.0*fabs(Kmax[i1])) && (-Kmin[i2] > 2.0*fabs(Kmax[i2]));
				//r->strong = (-r->k > fabs(r->k2));
				//r->k = fabs(r->k) - fabs(k2);
			}
		}
	}

	//connect ridge
	int* crossN = new int[faceN];
	bool* visit = new bool[faceN];
	for(i=0; i<faceN; i++){
		crossN[i] = 0;
		int j;  for(i = 0; j<3; j++){
			if(ridgeP[i][j] != NULL)
				crossN[i]++;
		}
		visit[i] = false;
	}

	//add special points
	/*
	for(i=0; i<faceN; i++){
		if(crossN[i] != 1)
			continue;

		if(face_link[i][0] < 0 || face_link[i][1] < 0 || face_link[i][2] < 0)
			continue;

		int j;
		if(ridgeP[i][0] != NULL)
			j = 0;
		else if(ridgeP[i][1] != NULL)
			j = 1;
		else
			j = 2;

		int ad = face_link[i][(j+1)%3];
		if(crossN[ad] != 0){
			int index;
			if(face_link[ad][0] == i)
				index = 0;
			else if(face_link[ad][1] == i)
				index = 1;
			else
				index = 2;

			RPoint* r = new RPoint;
			float* v1 = vertex[face[i][(j+1)%3]];
			float* v2 = vertex[face[i][(j+2)%3]];
			r->p[0] = 0.5f*(v1[0] + v2[0]);
			r->p[1] = 0.5f*(v1[1] + v2[1]);
			r->p[2] = 0.5f*(v1[2] + v2[2]);
			r->k = ridgeP[i][j]->k;
			r->strong = false;
			ridgeP[i][(j+1)%3] = ridgeP[ad][index] = r;
			
			crossN[ad]++;
			crossN[i]++;
		}

		ad = face_link[i][(j+2)%3];
		if(crossN[ad] != 0){
			int index;
			if(face_link[ad][0] == i)
				index = 0;
			else if(face_link[ad][1] == i)
				index = 1;
			else
				index = 2;

			RPoint* r = new RPoint;
			float* v1 = vertex[face[i][(j+2)%3]];
			float* v2 = vertex[face[i][j]];
			r->p[0] = 0.5f*(v1[0] + v2[0]);
			r->p[1] = 0.5f*(v1[1] + v2[1]);
			r->p[2] = 0.5f*(v1[2] + v2[2]);
			r->k = ridgeP[i][j]->k;
			r->strong = false;
			ridgeP[i][(j+2)%3] = ridgeP[ad][index] = r;
			
			crossN[ad]++;
			crossN[i]++;
		}
	}*/

	int* stack = new int[mesh->vertex_N];
	int top;
	mesh->ridge_S = new RHead;
	mesh->ridge_S->next = NULL;
	for(i=0; i<faceN; i++){
		if(crossN[i] < 2 || visit[i])
			continue;

		bool flag = false;
		int j;  for(i = 0; j< 3; j++){
			if(ridgeP[i][j] != NULL && ridgeP[i][j]->strong){
				flag = true;
				break;
			}
		}
		if(!flag)
			continue;

		visit[i] = true;
		RHead* head = new RHead;
		head->strength = 0;
		REdge* e = new REdge;
		e->next = NULL;

		stack[0] = i;
		top = 0;
		
		while(top >= 0){
			int current = stack[top--];

			if(crossN[current] == 1)
				continue;
			else if(crossN[current] == 2){
				RPoint *r1, *r2;
				r1 = NULL;
				int j;  for(i = 0; j<3; j++){
					if(ridgeP[current][j] != NULL){
						if(r1 == NULL){
							r1 = ridgeP[current][j];
							int ad = face_link[current][j];
							if(!visit[ad]){
								stack[++top] = ad;
								visit[ad] = true;
							}
						}
						else{
							r2 = ridgeP[current][j];
							int ad = face_link[current][j];
							if(!visit[ad]){
								stack[++top] = ad;
								visit[ad] = true;
							}
						}
					}
				}
				REdge* e1 = new REdge;
				e1->k = 0.5*(r1->k + r2->k);
				e1->p1[0] = r1->p[0];
				e1->p1[1] = r1->p[1];
				e1->p1[2] = r1->p[2];
				e1->p2[0] = r2->p[0];
				e1->p2[1] = r2->p[1];
				e1->p2[2] = r2->p[2];
				e1->next = e;
				e = e1;
				
				float vx = r1->p[0] - r2->p[0];
				float vy = r1->p[1] - r2->p[1];
				float vz = r1->p[2] - r2->p[2];
				double d = sqrt(vx*vx + vy*vy + vz*vz);
				head->strength += 0.5*(fabs(r1->k) + fabs(r2->k) 
					                   /*- fabs(r1->k2) - fabs(r2->k2)*/)*d;
			}
			else{
				RPoint** r = ridgeP[current];
				int j;  for(i = 0; j<3; j++){
					if(r[j] != NULL){
						int ad = face_link[current][j];
						if(!visit[ad]){
							stack[++top] = ad;
							visit[ad] = true;
						}
					}
				}

				REdge* e1 = new REdge;
				REdge* e2 = new REdge;
				REdge* e3 = new REdge;
				e1->next = e2;
				e2->next = e3;
				e3->next = e;
				e = e1;

				double k = (r[0]->k + r[1]->k + r[2]->k)/3.0;
				float cx = (r[0]->p[0] + r[1]->p[0] + r[2]->p[0])/3.0f;
				float cy = (r[0]->p[1] + r[1]->p[1] + r[2]->p[1])/3.0f;
				float cz = (r[0]->p[2] + r[1]->p[2] + r[2]->p[2])/3.0f;

				e1->k = k;
				e1->p1[0] = cx;
				e1->p1[1] = cy;
				e1->p1[2] = cz;
				e1->p2[0] = r[0]->p[0];
				e1->p2[1] = r[0]->p[1];
				e1->p2[2] = r[0]->p[2];
				
				e2->k = k;
				e2->p1[0] = cx;
				e2->p1[1] = cy;
				e2->p1[2] = cz;
				e2->p2[0] = r[1]->p[0];
				e2->p2[1] = r[1]->p[1];
				e2->p2[2] = r[1]->p[2];

				e3->k = k;
				e3->p1[0] = cx;
				e3->p1[1] = cy;
				e3->p1[2] = cz;
				e3->p2[0] = r[2]->p[0];
				e3->p2[1] = r[2]->p[1];
				e3->p2[2] = r[2]->p[2];
			}
		}

		head->head = e;
		head->next = mesh->ridge_S; 
		mesh->ridge_S = head;
	}
	for(i=0; i<faceN; i++){
		int j;  for(i = 0; j<3; j++){
			int ad = face_link[i][j];
			if(i > ad)
				continue;
			if(ridgeP[i][j] != NULL){
				delete ridgeP[i][j];
			}
		}
	}

	//connect ravine
	for(i=0; i<faceN; i++){
		crossN[i] = 0;
		int j;  for(i = 0; j<3; j++){
			if(ravineP[i][j] != NULL)
				crossN[i]++;
		}
		visit[i] = false;
	}

	//add special points
	/*
	for(i=0; i<faceN; i++){
		if(crossN[i] != 1)
			continue;

		if(face_link[i][0] < 0 || face_link[i][1] < 0 || face_link[i][2] < 0)
			continue;

		int j;
		if(ravineP[i][0] != NULL)
			j = 0;
		else if(ravineP[i][1] != NULL)
			j = 1;
		else
			j = 2;

		int ad = face_link[i][(j+1)%3];
		if(crossN[ad] == 1){
			int index;
			if(face_link[ad][0] == i)
				index = 0;
			else if(face_link[ad][1] == i)
				index = 1;
			else
				index = 2;

			RPoint* r = new RPoint;
			float* v1 = vertex[face[i][(j+1)%3]];
			float* v2 = vertex[face[i][(j+2)%3]];
			r->p[0] = 0.5f*(v1[0] + v2[0]);
			r->p[1] = 0.5f*(v1[1] + v2[1]);
			r->p[2] = 0.5f*(v1[2] + v2[2]);
			r->k = ravineP[i][j]->k;
			r->strong = false;
			ravineP[i][(j+1)%3] = ravineP[ad][index] = r;
			
			crossN[ad]++;
			crossN[i]++;
		}

		ad = face_link[i][(j+2)%3];
		if(crossN[ad] == 1){
			int index;
			if(face_link[ad][0] == i)
				index = 0;
			else if(face_link[ad][1] == i)
				index = 1;
			else
				index = 2;

			RPoint* r = new RPoint;
			float* v1 = vertex[face[i][(j+2)%3]];
			float* v2 = vertex[face[i][j]];
			r->p[0] = 0.5f*(v1[0] + v2[0]);
			r->p[1] = 0.5f*(v1[1] + v2[1]);
			r->p[2] = 0.5f*(v1[2] + v2[2]);
			r->k = ravineP[i][j]->k;
			r->strong = false;
			ravineP[i][(j+2)%3] = ravineP[ad][index] = r;
			
			crossN[ad]++;
			crossN[i]++;
		}
	}*/

	mesh->ravine_S = new RHead;
	mesh->ravine_S->next = NULL;
	for(i=0; i<faceN; i++){
		if(crossN[i] < 2 || visit[i])
			continue;

		
		bool flag = false;
		int j;  for(i = 0; j< 3; j++){
			if(ravineP[i][j] != NULL && ravineP[i][j]->strong){
				flag = true;
				break;
			}
		}
		if(!flag)
			continue;

		visit[i] = true;
		RHead* head = new RHead;
		head->strength = 0;
		REdge* e = new REdge;
		e->next = NULL;

		stack[0] = i;
		top = 0;
		
		while(top >= 0){
			int current = stack[top--];

			if(crossN[current] == 1)
				continue;
			else if(crossN[current] == 2){
				RPoint *r1, *r2;
				r1 = NULL;
				int j;  for(i = 0; j<3; j++){
					if(ravineP[current][j] != NULL){
						if(r1 == NULL){
							r1 = ravineP[current][j];
							int ad = face_link[current][j];
							if(!visit[ad]){
								stack[++top] = ad;
								visit[ad] = true;
							}
						}
						else{
							r2 = ravineP[current][j];
							int ad = face_link[current][j];
							if(!visit[ad]){
								stack[++top] = ad;
								visit[ad] = true;
							}
						}
					}
				}
				REdge* e1 = new REdge;
				e1->k = 0.5*(r1->k + r2->k);
				e1->p1[0] = r1->p[0];
				e1->p1[1] = r1->p[1];
				e1->p1[2] = r1->p[2];
				e1->p2[0] = r2->p[0];
				e1->p2[1] = r2->p[1];
				e1->p2[2] = r2->p[2];
				e1->next = e;
				e = e1;
				
				float vx = r1->p[0] - r2->p[0];
				float vy = r1->p[1] - r2->p[1];
				float vz = r1->p[2] - r2->p[2];
				double d = sqrt(vx*vx + vy*vy + vz*vz);
				head->strength += 0.5*(fabs(r1->k) + fabs(r2->k) 
					                   /*- fabs(r1->k2) - fabs(r2->k2)*/)*d;
			}
			else{
				RPoint** r = ravineP[current];
				int j;  for(i = 0; j<3; j++){
					if(r[j] != NULL){
						int ad = face_link[current][j];
						if(!visit[ad]){
							stack[++top] = ad;
							visit[ad] = true;
						}
					}
				}

				REdge* e1 = new REdge;
				REdge* e2 = new REdge;
				REdge* e3 = new REdge;
				e1->next = e2;
				e2->next = e3;
				e3->next = e;
				e = e1;

				double k = (r[0]->k + r[1]->k + r[2]->k)/3.0;
				float cx = (r[0]->p[0] + r[1]->p[0] + r[2]->p[0])/3.0f;
				float cy = (r[0]->p[1] + r[1]->p[1] + r[2]->p[1])/3.0f;
				float cz = (r[0]->p[2] + r[1]->p[2] + r[2]->p[2])/3.0f;

				e1->k = k;
				e1->p1[0] = cx;
				e1->p1[1] = cy;
				e1->p1[2] = cz;
				e1->p2[0] = r[0]->p[0];
				e1->p2[1] = r[0]->p[1];
				e1->p2[2] = r[0]->p[2];
				
				e2->k = k;
				e2->p1[0] = cx;
				e2->p1[1] = cy;
				e2->p1[2] = cz;
				e2->p2[0] = r[1]->p[0];
				e2->p2[1] = r[1]->p[1];
				e2->p2[2] = r[1]->p[2];

				e3->k = k;
				e3->p1[0] = cx;
				e3->p1[1] = cy;
				e3->p1[2] = cz;
				e3->p2[0] = r[2]->p[0];
				e3->p2[1] = r[2]->p[1];
				e3->p2[2] = r[2]->p[2];
			}
		}

		head->head = e;
		head->next = mesh->ravine_S; 
		mesh->ravine_S = head;
	}
	for(i=0; i<faceN; i++){
		int j;  for(i = 0; j<3; j++){
			int ad = face_link[i][j];
			if(i > ad)
				continue;
			if(ravineP[i][j] != NULL){
				delete ravineP[i][j];
			}
		}
	}

	delete[] visit;
	delete[] stack;
	delete[] crossN;
	delete[] ridgeP;
	delete[] ravineP;
}
