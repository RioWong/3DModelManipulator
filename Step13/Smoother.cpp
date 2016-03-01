/************************************************************************
Yutaka Ohtake 

Smoother.cpp
Smoothing methods

Copyright (c) 1999-2001 
The University of Aizu. All Rights Reserved.
************************************************************************/
/*
#include "stdafx.h"
#include "MeshEditor.h"

#include "FieldEnergy.h"
#include "EnergyMinimizer.h"
#include "Membrace.h"
#include "ThinPlate.h"
#include "NormalError.h"
#include "PBCG.h"
*/
#include "Smoother.h"


#define E 2.7182818284590
#define PI 3.14159265

//////////////////////////////////////////////////////////////////////
// \’z/Á–Å
//////////////////////////////////////////////////////////////////////

Smoother::Smoother()
{
	T = NULL;
	sigma = NULL;
}

Smoother::~Smoother()
{
	if(sigma != NULL)
		delete[] sigma;
	if(T != NULL)
		delete[] T;
}

//void Smoother::LaplacianFlow(int iter, float dt)
//{
//	/*
//	int n = mesh->vertex_N;
//
//	
//	Membrace e;
//	e.link = mesh->vertex_link_v;
//	e.degree = mesh->degree_v;
//	e.vertex = mesh->vertex;
//	e.n = mesh->vertex_N;
//	EnergyMinimizer mini;
//	mini.energy = &e;
//	
//	
//	ThinPlate e;
//	e.mesh = mesh;
//	EnergyMinimizer mini;
//	mini.energy = &e;
//	
//	int iter1;
//	float fret;
//	float *p = new float[3*n];
//	int i;  for(i = 0; i<n; i++){
//		p[3*i] = mesh->vertex[i][0];
//		p[3*i+1] = mesh->vertex[i][1];
//		p[3*i+2] = mesh->vertex[i][2];
//	}
//
//	mini.frprmn(p, 3*n, 0.000001f, &iter1, &fret);
//
//	for(i=0; i<n; i++){
//		mesh->vertex[i][0] = p[3*i];
//		mesh->vertex[i][1] = p[3*i+1];
//		mesh->vertex[i][2] = p[3*i+2];
//	}*/
//
//	
//	int vertex_N = mesh->vertex_N;
//	//L is Lplacian vectors
//	double (*L)[3] = new double[vertex_N][3];
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			//Boundary points are fixed
//			if(mesh->isBound[j]){
//				L[j][0] = L[j][1] = L[j][2] = 0;
//			}	
//			else{
//				mesh->laplacian(j, L[j]);
//			}
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			float *p = mesh->vertex[j];
//			//replace each vertex
//			p[0] += (float)(dt*L[j][0]);
//			p[1] += (float)(dt*L[j][1]);
//			p[2] += (float)(dt*L[j][2]);
//		}
//	}
//	delete[] L;
//	
//}

//void Smoother::meanCurvatureFlow(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	//Hn is mean curvature normals
//	double (*Hn)[3] = new double[vertex_N][3];
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			//Boundary points are fixed
//			if(mesh->isBound[j]){
//				Hn[j][0] = Hn[j][1] = Hn[j][2] = 0;
//			}
//			else{
//				//compute mean curvature normal
//				mesh->meanCurvatureNormal(j, Hn[j]);
//			}
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			float *p = mesh->vertex[j];
//			//replace each vertex
//			p[0] -= (float)(dt*Hn[j][0]);
//			p[1] -= (float)(dt*Hn[j][1]);
//			p[2] -= (float)(dt*Hn[j][2]);
//		}
//	}
//	delete[] Hn;
//}

void Smoother::setMeshData(MeshData *mesh)
{
	this->mesh = mesh;
}

//void Smoother::GaussianNoise(float size)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal)[3] = mesh->normal;
//
//	float noise;
//	srand(100);
//	int i;  for(i = 0; i<vertex_N; i++){
//		double r1, r2, r;
//		do{
//			r1 = 2.0*((double)rand())/RAND_MAX - 1.0;
//			r2 = 2.0*((double)rand())/RAND_MAX - 1.0;
//			r = r1*r1 + r2*r2;
//		}while(r >= 1.0 || r == 0.0);
//		double f = sqrt(-2.0*log(r)/r);
//		
//if(vertex[i][1] > 0)
//r1 = 0;
//
//		noise = (float)(size*r1*f);
//		vertex[i][0] += normal[i][0]*noise;
//		vertex[i][1] += normal[i][1]*noise;
//		vertex[i][2] += normal[i][2]*noise;
//
//		if(++i == vertex_N)
//			return;
//		else{
//			
//if(vertex[i][1] > 0)
//r2 = 0;
//
//			noise = (float)(size*r2*f);		
//			vertex[i][0] += normal[i][0]*noise;
//			vertex[i][1] += normal[i][1]*noise;
//			vertex[i][2] += normal[i][2]*noise;
//		}
//	}
//}

//void Smoother::computeNRingMeanCurvature(int n)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//
//	if(T != NULL)
//		delete[] T;
//
//	T = new float[vertex_N][3];
//
//	int i;  for(i = 0; i<vertex_N; i++){
//
//		int **nei, *nei_N;
//		mesh->getConsistentNRings(i, n, nei, nei_N);
//
//		float total_w = 0;
//		float total_Hn[3];
//		total_Hn[0] = total_Hn[1] = total_Hn[2] = 0;
//
//		for(int j=n-1; j<n; j++){
//
//			if(nei_N[j] == 0 || mesh->isBound[nei[j][0]])
//				continue;
//
//			float A_j = 0;
//			float Hn_j[3];
//			Hn_j[0] = Hn_j[1] = Hn_j[2] = 0;
//
//			for(int k=0; k<nei_N[j]; k+=j+1){
//				int s = nei[j][k];
//				int t = nei[j][(k+j+1)%nei_N[j]];
//
//				//vectors of edges of triangles
//				float PQ1[3], PQ2[3], QQ[3];
//				MeshData::VEC(PQ1, vertex[i], vertex[s]);
//				MeshData::VEC(PQ2, vertex[i], vertex[t]);
//				MeshData::VEC(QQ, vertex[s], vertex[t]);
//
//				//normal vector of triangle
//				float normal_t[3];
//				MeshData::CROSS(normal_t, PQ1, PQ2); 
//
//				//area of triangle
//				float A_k = (float)MeshData::AREA(vertex[i], vertex[s], vertex[t]);
//		
//				if(A_k != 0){
//					A_j += A_k;
//
//					float dot1 = (float)MeshData::DOT(PQ1, QQ);
//					float dot2 = -(float)MeshData::DOT(PQ2,QQ);
//
//					float cot1 = dot2/A_k;
//					float cot2 = dot1/A_k;
//
//					Hn_j[0] += cot1*PQ1[0] + cot2*PQ2[0];
//					Hn_j[1] += cot1*PQ1[1] + cot2*PQ2[1];
//					Hn_j[2] += cot1*PQ1[2] + cot2*PQ2[2];
//				}
//			}
//			if(A_j != 0){
//				Hn_j[0] /= A_j;
//				Hn_j[1] /= A_j;
//				Hn_j[2] /= A_j;
//
//				total_w += 1.0f; ///A_j;
//
//				total_Hn[0] += Hn_j[0]; ///A_j;
//				total_Hn[1] += Hn_j[1]; ///A_j;
//				total_Hn[2] += Hn_j[2]; ///A_j;
//			}
//		}
//	
//		for(j=0; j<n; j++)
//			delete[] nei[j];
//		delete[] nei;
//		delete[] nei_N;
//
//		if(total_w != 0){
//			T[i][0] = total_Hn[0]/total_w;
//			T[i][1] = total_Hn[1]/total_w;
//			T[i][2] = total_Hn[2]/total_w;
//		}
//		else
//			T[i][0] = T[i][1] = T[i][2] = 0;
//	}
//}

//void Smoother::meanCurvatureFlow12Ring(int iter, float dt, float e)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	//Hn is mean curvature normals
//	double (*Hn)[3] = new double[vertex_N][3];
//
//	// index of 2-ring
//	int **nei2, *nei2_N;
//	nei2 = (int **)new int[vertex_N];
//	nei2_N = new int[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++){
//		int **nei, *nei_N;
//		mesh->getConsistentNRings(i, 2, nei, nei_N);
//
//		nei2_N[i] = nei_N[1];
//		nei2[i] = nei[1];
//
//		delete[] nei[0];
//		delete[] nei;
//		delete[] nei_N;
//	}
//
//	for(i=0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			//Boundary points are fixed
//			if(mesh->isBound[j]){
//				Hn[j][0] = Hn[j][1] = Hn[j][2] = 0;
//			}
//			else{
//				//compute mean curvature normal
//				mesh->meanCurvatureNormal(j, Hn[j]);
//			
//				if(nei2_N[j] == 0 || mesh->isBound[nei2[j][0]])
//					continue;
//
//				float A_2 = 0;
//				float Hn_2[3];
//				Hn_2[0] = Hn_2[1] = Hn_2[2] = 0;
//
//				for(int k=0; k<nei2_N[j]; k += 2){
//					int s = nei2[j][k];
//					int t = nei2[j][(k+2)%nei2_N[j]];
//
//					//vectors of edges of triangles
//					float PQ1[3], PQ2[3], QQ[3];
//					MeshData::VEC(PQ1, vertex[j], vertex[s]);
//					MeshData::VEC(PQ2, vertex[j], vertex[t]);
//					MeshData::VEC(QQ, vertex[s], vertex[t]);
//
//					//normal vector of triangle
//					float normal_t[3];
//					MeshData::CROSS(normal_t, PQ1, PQ2); 
//
//					//area of triangle
//					float A_k = (float)MeshData::AREA(vertex[j], vertex[s], vertex[t]);
//		
//					if(A_k != 0){
//						A_2 += A_k;
//
//						float dot1 = (float)MeshData::DOT(PQ1, QQ);
//						float dot2 = -(float)MeshData::DOT(PQ2,QQ);
//
//						float cot1 = dot2/A_k;
//						float cot2 = dot1/A_k;
//
//						Hn_2[0] += cot1*PQ1[0] + cot2*PQ2[0];
//						Hn_2[1] += cot1*PQ1[1] + cot2*PQ2[1];
//						Hn_2[2] += cot1*PQ1[2] + cot2*PQ2[2];
//					}
//				}
//				if(A_2 != 0){
//					Hn_2[0] /= A_2;
//					Hn_2[1] /= A_2;
//					Hn_2[2] /= A_2;
//				}
//
//
//				double H2 = MeshData::DOT(Hn_2, Hn_2);
//				double H12 = MeshData::DOT(Hn[j], Hn_2);
//			
//				if(H12 < 0){
//				}
//				else if(H12 > (1.0f+e)*H2){
//					Hn[j][0] *= (1.0f-H2/H12);
//					Hn[j][1] *= (1.0f-H2/H12);
//					Hn[j][2] *= (1.0f-H2/H12);
//				}
//				//else if(H12 < (1.0f-e)*H2){
//					//Hn[j][0] *= (H2/H12-1.0f);
//					//Hn[j][1] *= (H2/H12-1.0f);
//					//Hn[j][2] *= (H2/H12-1.0f);
//				//}
//				else{
//					Hn[j][0] = Hn[j][1] = Hn[j][2] = 0;
//				}
//			}
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			float *p = vertex[j];
//			//replace each vertex
//			p[0] -= (float)(dt*Hn[j][0]);
//			p[1] -= (float)(dt*Hn[j][1]);
//			p[2] -= (float)(dt*Hn[j][2]);
//		}
//	}
//	delete[] Hn;
//
//	for(i=0; i<vertex_N; i++){
//		delete[] nei2[i];
//	}
//	delete[] nei2;
//	delete[] nei2_N;
//}

//void Smoother::LaplacianWithRidge(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	BOOL *isFeature = new BOOL[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++)
//		isFeature[i] = (mesh->isRidge[i] || mesh->isRavine[i]);
//
//	//L is Lplacian vectors
//	double (*L)[3] = new double[vertex_N][3];
//
//	for(i=0; i<iter; i++){
//		mesh->computeFaceNormal();
//		mesh->computeNormal();
//		int j;  for(i = 0; j<vertex_N; j++){
//			L[j][0] = L[j][1] = L[j][2] = 0;
//			if(!mesh->isBound[j]){
//				int *nei, nei_N;
//				mesh->getConsistent1Ring(j, nei, nei_N);
//				if(nei_N == 0){
//					L[j][0] = L[j][1] = L[j][2] = 0;
//					continue;
//				}
//
//				int cnt = 0;
//				if(isFeature[j]){
//					for(int k=0; k<nei_N; k++){
//						if(isFeature[nei[k]]){
//							L[j][0] += vertex[nei[k]][0];
//							L[j][1] += vertex[nei[k]][1];
//							L[j][2] += vertex[nei[k]][2];
//							cnt++;
//						}
//					}
//				}
//
//				if(cnt != 0){
//					L[j][0] = L[j][0]/cnt - vertex[j][0];
//					L[j][1] = L[j][1]/cnt - vertex[j][1];
//					L[j][2] = L[j][2]/cnt - vertex[j][2];
//				}
//				else
//					mesh->laplacian(j, L[j]);
//
//				if(cnt == 1){
//					mesh->laplacian(j, L[j]);
//					double inner = MeshData::DOT(mesh->normal[j], L[j]);
//					L[i][0] = inner*mesh->normal[j][0];
//					L[i][1] = inner*mesh->normal[j][1];
//					L[i][2] = inner*mesh->normal[j][2];
//				}
//			}
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			float *p = mesh->vertex[j];
//			//replace each vertex
//			p[0] += (float)(dt*L[j][0]);
//			p[1] += (float)(dt*L[j][1]);
//			p[2] += (float)(dt*L[j][2]);
//		}
//	}
//	delete[] L;
//	delete[] isFeature;
//}

//void Smoother::smoothTmaxTmin(int iter, float T)
//{
//	int vertex_N = mesh->vertex_N;
//	double (*t_max)[3] = mesh->t_max;
//	double (*t_min)[3] = mesh->t_min;
//	float (*vertex)[3] = mesh->vertex;
//
//	double (*t_max1)[3] = new double[vertex_N][3];
//	double (*t_min1)[3] = new double[vertex_N][3];
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			int *nei, nei_N;
//			mesh->getConsistent1Ring(j, nei, nei_N);
//			if(nei_N == 0){
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//				continue;
//			}
//			
//			t_max1[j][0] = t_max1[j][1] = t_max1[j][2] = 0;
//			t_min1[j][0] = t_min1[j][1] = t_min1[j][2] = 0;
//			for(int k=0; k<nei_N; k++){
//				int pair = nei[k];
//				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
//							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
//							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
//
//				double dot = MeshData::DOT(t_max[j], t_max[pair]);
//				double dk = acos(fabs(dot));
//				double w = exp(-dk*dk/d2/(T*T));
//				if(dot < 0)
//					w = -w;
//				t_max1[j][0] += w*t_max[pair][0];
//				t_max1[j][1] += w*t_max[pair][1];
//				t_max1[j][2] += w*t_max[pair][2];
//				
//				dot = MeshData::DOT(t_min[j], t_min[pair]);
//				dk = acos(fabs(dot));
//				w = exp(-dk*dk/d2/(T*T));
//				if(dot < 0)
//					w = -w;
//				t_min1[j][0] += w*t_min[pair][0];
//				t_min1[j][1] += w*t_min[pair][1];
//				t_min1[j][2] += w*t_min[pair][2];
//			}
//
//			double len = MeshData::LENGTH(t_max1[j]);
//			if((float)len != 0){
//				t_max1[j][0] /= len;
//				t_max1[j][1] /= len;
//				t_max1[j][2] /= len;
//			}
//			else{
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//			}
//
//			len = MeshData::LENGTH(t_min1[j]);
//			if((float)len != 0){
//				t_min1[j][0] /= len;
//				t_min1[j][1] /= len;
//				t_min1[j][2] /= len;
//			}
//			else{
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//			}
//		}
//		for(j=0; j<vertex_N; j++){
//			t_max[j][0] = t_max1[j][0];
//			t_max[j][1] = t_max1[j][1];
//			t_max[j][2] = t_max1[j][2];
//			
//			t_min[j][0] = t_min1[j][0];
//			t_min[j][1] = t_min1[j][1];
//			t_min[j][2] = t_min1[j][2];
//		}
//	}
//
//	delete[] t_max1;
//	delete[] t_min1;
//
//	float (*normal)[3] = mesh->normal;
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, normal[i], t_max[i]);
//		if(MeshData::DOT(n, t_min[i])< 0){
//			t_min[i][0] = -t_min[i][0];
//			t_min[i][1] = -t_min[i][1];
//			t_min[i][2] = -t_min[i][2];
//		}	
//	}
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, t_max[i], t_min[i]);
//		double len = MeshData::LENGTH(n);
//		if((float)len != 0){
//			normal[i][0] = (float)(n[0]/len);
//			normal[i][1] = (float)(n[1]/len);
//			normal[i][2] = (float)(n[2]/len);
//		}
//		if(fabs(mesh->k_max[i]) > fabs(mesh->k_min[i])){
//			MeshData::CROSS(n, normal[i], t_max[i]);
//			double len = MeshData::LENGTH(n);
//			if((float)len != 0){
//				t_min[i][0] = n[0]/len;
//				t_min[i][1] = n[1]/len;
//				t_min[i][2] = n[2]/len;
//			}
//		}
//		else{
//			MeshData::CROSS(n, t_min[i], normal[i]);
//			double len = MeshData::LENGTH(n);
//			if((float)len != 0){
//				t_max[i][0] = n[0]/len;
//				t_max[i][1] = n[1]/len;
//				t_max[i][2] = n[2]/len;
//			}
//		}
//	}
//}

//void Smoother::smoothKmaxKmin(int iter, float T)
//{
//	int vertex_N = mesh->vertex_N;
//	double *k_max = mesh->k_max;
//	double *k_min = mesh->k_min;
//	float (*vertex)[3] = mesh->vertex;
//
//	double *k_max1 = new double[vertex_N];
//	double *k_min1 = new double[vertex_N];
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			int *nei, nei_N;
//			mesh->getConsistent1Ring(j, nei, nei_N);
//			if(nei_N == 0){
//				k_max1[j] = k_max[j];
//				k_min1[j] = k_min[j];
//				continue;
//			}
//			
//			double total_max = 0;
//			double total_min = 0;
//			k_max1[j] = k_min1[j] = 0;
//			for(int k=0; k<nei_N; k++){
//				int pair = nei[k];
//				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
//							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
//							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
//
//				double dk = k_max[j] - k_max[pair];
//				double w = exp(-dk*dk/d2/(T*T));
//				//float w = exp(-d2/(T*T));
//				total_max += w;
//				k_max1[j] += w*k_max[pair];
//
//				dk = k_min[j] - k_min[pair];
//				w = exp(-dk*dk/d2/(T*T));
//				//w = exp(-d2/(T*T));
//				total_min += w;
//				k_min1[j] += w*k_min[pair];
//			}
//			if(total_max != 0)
//				k_max1[j] /= total_max;
//			else
//				k_max1[j] = k_max[j];
//
//			if(total_min != 0)
//				k_min1[j] /= total_min;
//			else
//				k_min1[j] = k_min[j];
//		}
//		for(j=0; j<vertex_N; j++){
//			k_max[j] = k_max1[j];
//			k_min[j] = k_min1[j];
//		}
//	}
//
//	delete[] k_max1;
//	delete[] k_min1;
//}

//void Smoother::smoothTagEdges(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal)[3];
//	int *dist = new int[vertex_N];
//	int max_dist = 3;
//
//	Node** ridge = mesh->ridge_edge;
//	Node** ravine = mesh->ravine_edge;
//	Node** tag = new Node*[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++){
//		tag[i] = new Node;
//		for(Node* current = ridge[i]; current->next!=NULL; current=current->next)
//			tag[i]->append(current->v, -1);
//		for(Node * current = ravine[i]; current->next!=NULL; current=current->next)
//			tag[i]->append(current->v, -1);
//		if(tag[i]->next!=NULL)
//			dist[i] = 0;
//		else
//			dist[i] = 10000;
//	}
//
//	for(i=0; i<max_dist; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(dist[j] == i){
//				int *nei, nei_N;
//				mesh->get1Ring(j, nei, nei_N);
//				if(nei_N == 0)
//					continue;
//				for(int k=0; k<nei_N; k++)
//					dist[nei[k]] = min(dist[nei[k]], i+1);
//			}
//		}
//	}
//
//
//	double (*L)[3] = new double[vertex_N][3];
//	double (*BL)[3] = new double[vertex_N][3];
//	for(i=0; i<iter; i++){
//		//smooth tag curve 
//		mesh->computeFaceNormal();
//		mesh->computeNormal();
//		normal = mesh->normal;
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(dist[j] != 0)
//				continue;
//			L[j][0] = L[j][1] = L[j][2] = 0;
//			int m = 0;
//			for(Node* current=tag[j]; current->next!=NULL; current=current->next){
//				m++;
//				L[j][0] += vertex[current->v][0];
//				L[j][1] += vertex[current->v][1];
//				L[j][2] += vertex[current->v][2];
//			}
//			if(m == 2){
//				L[j][0] = L[j][0]/(float)m - vertex[j][0];
//				L[j][1] = L[j][1]/(float)m - vertex[j][1];
//				L[j][2] = L[j][2]/(float)m - vertex[j][2];
//			
//				double inner = MeshData::DOT(L[j], normal[j]);
//				L[j][0] = L[j][0] - inner*normal[j][0];
//				L[j][1] = L[j][1] - inner*normal[j][1];
//				L[j][2] = L[j][2] - inner*normal[j][2];
//			}
//			else{
//				L[j][0] = L[j][1] = L[j][2] = 0;
//			}
//		}
//		for(j=0; j<vertex_N; j++){
//			if(dist[j] == 0){
//				vertex[j][0] += (float)(dt*L[j][0]);
//				vertex[j][1] += (float)(dt*L[j][1]);
//				vertex[j][2] += (float)(dt*L[j][2]);
//			}
//		}
//
//		for(j=0; j<vertex_N; j++){
//			L[j][0] = L[j][1] = L[j][2] = 0;
//			if(dist[j] > max_dist)
//				continue;
//
//			int *nei, nei_N;
//			mesh->get1Ring(j, nei, nei_N);
//			if(nei_N == 0)
//				continue;
//
//			for(int k=0; k<nei_N; k++){
//				L[j][0] += vertex[nei[k]][0];
//				L[j][1] += vertex[nei[k]][1];
//				L[j][2] += vertex[nei[k]][2];
//			}
//		
//			L[j][0] = L[j][0]/(float)nei_N - vertex[j][0];
//			L[j][1] = L[j][1]/(float)nei_N - vertex[j][1];
//			L[j][2] = L[j][2]/(float)nei_N - vertex[j][2];
//		}
//
//		for(j=0; j<vertex_N; j++){
//			BL[j][0] = BL[j][1] = BL[j][2] = 0;
//			if(dist[j] == 0 || dist[j] >= max_dist)
//				continue;
//	
//			int *nei, nei_N;
//			mesh->get1Ring(j, nei, nei_N);
//			if(nei_N == 0)
//				continue;
//
//			for(int k=0; k<nei_N; k++){
//				BL[j][0] += L[nei[k]][0];
//				BL[j][1] += L[nei[k]][1];
//				BL[j][2] += L[nei[k]][2];
//			}
//		
//			BL[j][0] = BL[j][0]/(float)nei_N - L[j][0];
//			BL[j][1] = BL[j][1]/(float)nei_N - L[j][1];
//			BL[j][2] = BL[j][2]/(float)nei_N - L[j][2];
//			/*
//			int m = 0;
//			for(Node* current=tag[j]; current->next!=NULL; current=current->next){
//				m++;
//				BL[j][0] += L[current->v][0];
//				BL[j][1] += L[current->v][1];
//				BL[j][2] += L[current->v][2];
//			}
//			if(m!=0){
//				BL[j][0] = BL[j][0]/(float)m - L[j][0];
//				BL[j][1] = BL[j][1]/(float)m - L[j][1];
//				BL[j][2] = BL[j][2]/(float)m - L[j][2];
//			}
//			*/
//		}
//
//		for(j=0; j<vertex_N; j++){
//			if(!mesh->isBound[j]){
//				vertex[j][0] -= (float)(dt*BL[j][0]);
//				vertex[j][1] -= (float)(dt*BL[j][1]);
//				vertex[j][2] -= (float)(dt*BL[j][2]);
//			}
//		}
//	}
//	delete[] L;
//	delete[] BL;
//	delete[] dist;
//	for(i=0; i<vertex_N; i++)
//		delete tag[i];
//	delete tag;
//}

//void Smoother::smoothNormal(float T, int f_type, int d_type)
//{
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal)[3] = mesh->normal_f;
//	int (*face)[3] = mesh->face;
//	int face_N = mesh->face_N;
//	float (*tmp_normal)[3] = new float[face_N][3];
//	int (*ad_face)[3] = mesh->face_link_E;
//	int i;  for(i = 0; i<face_N; i++){
//		double total_w = 0;
//		double sum_n[3];
//		sum_n[0] = sum_n[1] = sum_n[2] = 0;
//		double c[3];
//		c[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3.0;
//		c[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3.0;
//		c[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3.0;
//		int j;  for(i = 0; j<3; j++){
//			int pair = ad_face[i][j];
//			if(pair < 0)
//				continue;
//		
//			double c1[3];
//			c1[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3.0;
//			c1[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3.0;
//			c1[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3.0;
//			double d = MeshData::DIST(c,c1);
//
//			double area = MeshData::AREA(vertex[face[pair][0]], vertex[face[pair][1]], vertex[face[pair][2]]);
//
//			double dot = MeshData::DOT(normal[i], normal[pair]);
//			if(dot > 1)
//				dot = 1;
//			else if(dot < -1)
//				dot = -1;
//
//			double k = acos(dot);
//			if(d_type == 0)
//				k /= d;
//
//			double w;
//			if(f_type == 1)
//				w = area*exp(-T*k*k);
//			else if(f_type == 0)
//				w = area/(1+T*k*k);
//			else 
//				w = area;
//				
//			total_w += w;
//			sum_n[0] += w*normal[pair][0];
//			sum_n[1] += w*normal[pair][1];
//			sum_n[2] += w*normal[pair][2];
//		}
//		if((float)total_w != 0){
//			tmp_normal[i][0] = (float)(sum_n[0]/total_w);
//			tmp_normal[i][1] = (float)(sum_n[1]/total_w);
//			tmp_normal[i][2] = (float)(sum_n[2]/total_w);
//
//			double len = MeshData::LENGTH(tmp_normal[i]);
//			if((float)len != 0){
//				tmp_normal[i][0] = (float)(tmp_normal[i][0]/len);
//				tmp_normal[i][1] = (float)(tmp_normal[i][1]/len);
//				tmp_normal[i][2] = (float)(tmp_normal[i][2]/len);
//			}
//			else{
//				tmp_normal[i][0] = normal[i][0];
//				tmp_normal[i][1] = normal[i][1];
//				tmp_normal[i][2] = normal[i][2];
//			}
//		}
//		else{
//			tmp_normal[i][0] = normal[i][0];
//			tmp_normal[i][1] = normal[i][1];
//			tmp_normal[i][2] = normal[i][2];
//		}
//	}
//	delete[] normal;
//	mesh->normal_f = tmp_normal;
//}

//void Smoother::moveNormal(float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	int (*face)[3] = mesh->face;
//	float (*normal_f)[3] = mesh->normal_f;
//	double (*dp)[3] = new double[vertex_N][3];
//	int i;  for(i = 0; i<vertex_N; i++){
//		if(mesh->isBound[i]){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//
//		int *nei, *nei_v, nei_N;
//		mesh->getConsistent1Ring(i, nei_v, nei, nei_N);
//		if(nei_N == 0){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//
//		double total_w = 0;
//
//		double pro[3];
//		pro[0] = pro[1] = pro[2] = 0;
//
//		BOOL flip_flag = false;
//
//		for(int j = 0; j<nei_N; j++){
//			int f = nei[j];
//
//			int i0 = face[f][0];
//			int i1 = face[f][1];
//			int i2 = face[f][2];
//
//			double o[3];
//			o[0] = (vertex[i0][0] + vertex[i1][0] + vertex[i2][0])/3.0;
//			o[1] = (vertex[i0][1] + vertex[i1][1] + vertex[i2][1])/3.0;
//			o[2] = (vertex[i0][2] + vertex[i1][2] + vertex[i2][2])/3.0;
//
//			double d[3];
//			MeshData::VEC(d, o, vertex[i]);
//			double inner = MeshData::DOT(normal_f[f], d);
//
//			double w = MeshData::AREA(vertex[i0], vertex[i1], vertex[i2]);
//			total_w += w;
//
//			pro[0] -= w*inner*normal_f[f][0];
//			pro[1] -= w*inner*normal_f[f][1];
//			pro[2] -= w*inner*normal_f[f][2];
//
//			float n[3];
//			float v1[3];
//			float v2[3];
//			MeshData::VEC(v1, vertex[i0], vertex[i1]);
//			MeshData::VEC(v2, vertex[i0], vertex[i2]);
//			MeshData::CROSS(n, v1, v2);
//			//if(MeshData::DOT(normal_f[f], n) > 0.7){
//				//flip_flag = true;
//				//break;
//				//TRACE("FLIP\n");
//			//}
//		}
//		if(flip_flag){
//			mesh->laplacian(i, dp[i]);
//			/*
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			int j;  for(i = 0; j<nei_N; j++){
//				dp[i][0] += vertex[nei_v[j]][0] - vertex[i][0];
//				dp[i][1] += vertex[nei_v[j]][1] - vertex[i][0];
//				dp[i][2] += vertex[nei_v[j]][2] - vertex[i][0];
//			}
//			dp[i][0] /= (float)nei_N;
//			dp[i][1] /= (float)nei_N;
//			dp[i][2] /= (float)nei_N;*/
//		}
//		else if(total_w != 0){
//			dp[i][0] = pro[0]/total_w;
//			dp[i][1] = pro[1]/total_w;
//			dp[i][2] = pro[2]/total_w;
//		}
//		else
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//	}
//
//	//mesh->computeFaceNormal();
//	//mesh->computeNormal();
//	for(i=0; i<vertex_N; i++){
//		vertex[i][0] += (float)(dt*dp[i][0]);
//		vertex[i][1] += (float)(dt*dp[i][1]);
//		vertex[i][2] += (float)(dt*dp[i][2]);
//	}
//	delete[] dp;
//}

//void Smoother::adaptiveSmoothing(int iter, float T, float rate, int int_type, int inner_iter, float int_step, int d_type, int ave_iter, int f_type, int L_deg)
//{
//	int i,j;
//	
//	if(mesh->face_link_V == NULL)
//		mesh->generateFaceLinkV();
//
//	for(i=0; i<iter; i++){
//		mesh->computeFaceNormal();
//		
//		if(rate != 0){
//			//float *sigma = new float[mesh->face_N];
//			//this->decideSigma(sigma, T, ave_iter);
//			//this->smoothNormalGaussian(ave_iter, sigma);
//			//delete[] sigma;
//			for(j=0; j<ave_iter; j++)
//				this->smoothNormalV(T,rate, f_type, d_type);
//			//float e = mesh->averageOfEdgeLength();
//			//computeOptimalGaussian(0.05, 1.0, 0.05);
//			//computeOptimalGaussian(0.25f*e, 2.5f*e, 0.25f*e, T*e*e);
//		}
//		else
//			for(j=0; j<ave_iter; j++)
//				this->smoothNormal(T, f_type, d_type);
//					
//		if(int_type == 0){
//			for(j=0; j<inner_iter; j++)
//				this->moveNormal(int_step);
//		}
//		else{
//			if(L_deg != 1){
//				for(j=0; j<inner_iter; j++)
//					this->moveNormalAreaD(int_step, L_deg);
//			}
//			else{
//				//for(j=0; j<inner_iter; j++)
//					//this->minimizeNormalErr();
//				//this->IntegrateNormalImplicit(int_step);
//				for(j=0; j<inner_iter; j++)
//					this->moveNormalAreaD(int_step);
//			}		
//		}
//	}
//}

//void Smoother::curvatureAlongMedian(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	double (*dp)[3] = new double[vertex_N][3];
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(mesh->isBound[j]){
//				dp[j][0] = dp[j][1] = dp[j][2] = 0;
//				continue;
//			}
//
//			double L[3], Hn[3];
//			mesh->laplacian(j, L);
//			mesh->meanCurvatureNormal(j, Hn);
//			Hn[0] *= -1;
//			Hn[1] *= -1;
//			Hn[2] *= -1;
//			double dot = L[0]*Hn[0] + L[1]*Hn[1] + L[2]*Hn[2];
//			double H2 = Hn[0]*Hn[0] + Hn[1]*Hn[1] + Hn[2]*Hn[2];
//			double L2 = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
//	
//			if(dot == 0){
//				dp[j][0] = L[0];
//				dp[j][1] = L[1];
//				dp[j][2] = L[2];
//			}
//			else if(dot > 0){//.1*sqrt(H2)*sqrt(L2)){
//				dp[j][0] = dt*(float)(H2*L[0]/dot);
//				dp[j][1] = dt*(float)(H2*L[1]/dot);
//				dp[j][2] = dt*(float)(H2*L[2]/dot);
//
//				
//				double lenL = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
//				double lenDp = dp[j][0]*dp[j][0] + dp[j][1]*dp[j][1] + dp[j][2]*dp[j][2];
//				if(lenL < lenDp){
//					dp[j][0] = L[0];
//					dp[j][1] = L[1];
//					dp[j][2] = L[2];
//				}
//				
//			}
//			else if(dot < 0){ //-0.1*sqrt(H2)*sqrt(L2)){
//				dp[j][0] = dt*(float)(2.0*Hn[0] - H2*L[0]/dot);
//				dp[j][1] = dt*(float)(2.0*Hn[1] - H2*L[1]/dot);
//				dp[j][2] = dt*(float)(2.0*Hn[2] - H2*L[2]/dot);
//				
//				double lenL = L[0]*L[0] + L[1]*L[1] + L[2]*L[2];
//				double lenDp = dp[j][0]*dp[j][0] + dp[j][1]*dp[j][1] + dp[j][2]*dp[j][2];
//				if(lenL < lenDp){
//					dp[j][0] = L[0] - 2.0*dot*Hn[0]/H2;
//					dp[j][1] = L[1] - 2.0*dot*Hn[1]/H2;
//					dp[j][2] = L[2] - 2.0*dot*Hn[2]/H2;
//				}
//			}
//			else{
//				dp[j][0] = dp[j][1] = dp[j][2] = 0;
//			}
//		}
//
//		for(j=0; j<vertex_N; j++){
//			vertex[j][0] += (float)dp[j][0];
//			vertex[j][1] += (float)dp[j][1];
//			vertex[j][2] += (float)dp[j][2];
//		}
//	}
//	delete[] dp;
//}

void Smoother::smoothNormalV(float T, float rate, int f_type, int d_type)
{
	float (*vertex)[3] = mesh->vertex;
	float (*normal)[3] = mesh->normal_f;
	int (*face)[3] = mesh->face;
	int face_N = mesh->face_N;
	float (*tmp_normal)[3] = new float[face_N][3];
	int (*ad_face)[3] = mesh->face_link_E;
	int **ad_face2 = mesh->face_link_V;
	int *ad_face2_N = mesh->face_link_V_N;

	int i;  for(i = 0; i<face_N; i++){
		double total_w = 0;
		double sum_n[3];
		sum_n[0] = sum_n[1] = sum_n[2] = 0;
		double c[3];
		c[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3.0;
		c[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3.0;
		c[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3.0;
		int j;  for(i = 0; j<3; j++){
			int pair = ad_face[i][j];
			if(pair < 0)
				continue;
		
			double c1[3];
			c1[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3.0;
			c1[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3.0;
			c1[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3.0;
			double d = MeshData::DIST(c,c1);
			if((float)d == 0)
				continue;

			double area = MeshData::AREA(vertex[face[pair][0]], vertex[face[pair][1]], vertex[face[pair][2]]);

			double dot = MeshData::DOT(normal[i], normal[pair]);
			if(dot > 1)
				dot = 1;
			else if(dot < -1)
				dot = -1;

			double k = acos(dot);
			if(d_type == 0)
				k /= d;

			double w;
			if(f_type == 1)
				w = area*exp(-T*k*k);
			else if(f_type == 0)
				w = area/(1.0+T*k*k);
			else 
				w = 1; //area;
				
			total_w += w;
			sum_n[0] += w*normal[pair][0];
			sum_n[1] += w*normal[pair][1];
			sum_n[2] += w*normal[pair][2];
		}
		for(j=0; j<ad_face2_N[i]; j++){
			int pair = ad_face2[i][j];
			if(pair < 0 || pair >= face_N)
				continue;
		
			double c1[3];
			c1[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3.0;
			c1[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3.0;
			c1[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3.0;
			double d = MeshData::DIST(c,c1);
			if((float)d == 0)
				continue;

			double area = MeshData::AREA(vertex[face[pair][0]], vertex[face[pair][1]], vertex[face[pair][2]]);
			
			double dot = MeshData::DOT(normal[i], normal[pair]);
			if(dot > 1)
				dot = 1;
			else if(dot < -1)
				dot = -1;

			double k = acos(dot);
			if(d_type == 0)
				k /= d;

			double w;
			if(f_type == 1)
				w = rate*area*exp(-T*k*k);
			else if(f_type == 0)
				w = rate*area/(1.0+T*k*k);
			else 
				w = rate; //rate*area;

			total_w += w;
			sum_n[0] += w*normal[pair][0];
			sum_n[1] += w*normal[pair][1];
			sum_n[2] += w*normal[pair][2];
		}
		if((float)total_w != 0){
			tmp_normal[i][0] = (float)(sum_n[0]/total_w);
			tmp_normal[i][1] = (float)(sum_n[1]/total_w);
			tmp_normal[i][2] = (float)(sum_n[2]/total_w);

			double len = MeshData::LENGTH(tmp_normal[i]);
			if((float)len != 0){
				tmp_normal[i][0] = (float)(tmp_normal[i][0]/len);
				tmp_normal[i][1] = (float)(tmp_normal[i][1]/len);
				tmp_normal[i][2] = (float)(tmp_normal[i][2]/len);
			}
			else{
				tmp_normal[i][0] = normal[i][0];
				tmp_normal[i][1] = normal[i][1];
				tmp_normal[i][2] = normal[i][2];
			}
		}
		else{
			tmp_normal[i][0] = normal[i][0];
			tmp_normal[i][1] = normal[i][1];
			tmp_normal[i][2] = normal[i][2];
		}
	}
	delete[] normal;
	mesh->normal_f = tmp_normal;
}

void Smoother::smoothNormalV(int iter, float T, float rate)
{
	/*
	float (*vertex)[3] = mesh->vertex;
	float (*normal)[3] = mesh->normal_f;
	int (*face)[3] = mesh->face;
	int face_N = mesh->face_N;
	float (*tmp_normal)[3] = new float[face_N][3];
	int (*ad_face)[3] = mesh->face_link_E;
	int **ad_face2 = mesh->face_link_V;
	int *ad_face2_N = mesh->face_link_V_N;

	for(int m=0; m<iter; m++){
		int i;  for(i = 0; i<face_N; i++){
			float total_w = 0;
			tmp_normal[i][0] = 0;
			tmp_normal[i][1] = 0;
			tmp_normal[i][2] = 0;
			double c[3];
			c[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3.0;
			c[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3.0;
			c[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3.0;
			int j;  for(i = 0; j<3; j++){
				int pair = ad_face[i][j];
				if(pair < 0)
					continue;
		
				double c1[3];
				c1[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3.0;
				c1[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3.0;
				c1[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3.0;
				double d = MeshData::DIST(c,c1);
				if((float)d == 0)
					continue;

				float area = MeshData::AREA(vertex[face[pair][0]], vertex[face[pair][1]], vertex[face[pair][2]]);

				float k = acos(MeshData::DOT(normal[i], normal[pair]))/d;
				float w = area*exp(-(k/T)*(k/T));
				total_w += w;
				tmp_normal[i][0] += w*normal[pair][0];
				tmp_normal[i][1] += w*normal[pair][1];
				tmp_normal[i][2] += w*normal[pair][2];
			}
			for(j=0; j<ad_face2_N[i]; j++){
				int pair = ad_face2[i][j];
				if(pair < 0 || pair >= face_N)
					continue;
		
				double c1[3];
				c1[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3.0f;
				c1[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3.0f;
				c1[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3.0f;
				double d = MeshData::DIST(c,c1);
				if((float)d == 0)
					continue;

			float area = MeshData::AREA(vertex[face[pair][0]], vertex[face[pair][1]], vertex[face[pair][2]]);

			float k = acos(MeshData::DOT(normal[i], normal[pair]))/d;
			float w = rate*area*exp(-(k/T)*(k/T));
			total_w += w;
			tmp_normal[i][0] += w*normal[pair][0];
			tmp_normal[i][1] += w*normal[pair][1];
			tmp_normal[i][2] += w*normal[pair][2];
			}
			if(total_w != 0){
				tmp_normal[i][0] /= total_w;
				tmp_normal[i][1] /= total_w;
				tmp_normal[i][2] /= total_w;

				float len = MeshData::LENGTH(tmp_normal[i]);
				if(len != 0){
					tmp_normal[i][0] /= len;
					tmp_normal[i][1] /= len;
					tmp_normal[i][2] /= len;
				}
				else{
					tmp_normal[i][0] = normal[i][0];
					tmp_normal[i][1] = normal[i][1];
					tmp_normal[i][2] = normal[i][2];
				}
			}
			else{
				tmp_normal[i][0] = normal[i][0];
				tmp_normal[i][1] = normal[i][1];
				tmp_normal[i][2] = normal[i][2];
			}
		}
	}
	int i;  for(i = 0; i<face_N; i++){
		double len = MeshData::LENGTH(tmp_normal[i]);
		if(len != 0){
			tmp_normal[i][0] /= len;
			tmp_normal[i][1] /= len;
			tmp_normal[i][2] /= len;
		}
		else{
			tmp_normal[i][0] = normal[i][0];
			tmp_normal[i][1] = normal[i][1];
			tmp_normal[i][2] = normal[i][2];
		}
	}
	delete[] normal;
	mesh->normal_f = tmp_normal;
	*/
}

//void Smoother::moveMinPosition()
//{
//	int vertex_N = mesh->vertex_N;
//	int face_N = mesh->face_N;
//	float (*vertex)[3] = mesh->vertex;
//	int (*face)[3] = mesh->face;
//	float (*normal)[3] = mesh->normal_f;
//	int *type = mesh->isBound;
//	int **link = mesh->vertex_link_f;
//	int *degree = mesh->degree_f;
//	double (*center)[3] = new double[face_N][3];
//
//	int i;  for(i = 0; i<face_N; i++)
//		mesh->faceCenter(center[i], i);
//
//	for(i=0; i<vertex_N; i++){
//		if(type[i])
//			continue;
//
//		double S[6] = {0,0,0,0,0};
//		double b[3] = {0,0,0};
//		int f = degree[i];
//		int* l = link[i];
//		double total_a = 0;
//		int j;  for(i = 0; j<f; j++){
//			double A[6];
//			double n[3];
//			n[0] = (double)normal[l[j]][0];
//			n[1] = (double)normal[l[j]][1];
//			n[2] = (double)normal[l[j]][2];
//			MATRIX(A, n);
//			double a = mesh->faceArea(l[j]);
//			total_a += a;
//			MAT_TIMES(A, a);
//			MAT_SUM(S, A);
//			double d = MeshData::DOT(n, center[l[j]]);
//			b[0] += a*d*n[0];
//			b[1] += a*d*n[1];
//			b[2] += a*d*n[2];
//		}
//		if((float)total_a !=0){
//			for(int k=0; k<6; k++)
//				S[k] /= total_a;
//			b[0] /= total_a;
//			b[1] /= total_a;
//			b[2] /= total_a;
//		}
//		double InS[6];
//		if(INVERSE(InS, S)){
//			double p[3];
//			MAT_BY_VEC(p, InS, b);
//			//zone
//			BOOL flag = true;
//			for(j=0; j<f; j++){
//				int i1 = l[j];
//				int i2 = l[(j+1)%f];
//				double v[3];
//				MeshData::VEC(v, center[i1], center[i2]);
//				double nm[3];
//				nm[0] = normal[i1][0] + normal[i2][0];
//				nm[1] = normal[i1][1] + normal[i2][1];
//				nm[2] = normal[i1][2] + normal[i2][2];
//				double n[3];
//				MeshData::CROSS(n, v, nm);
//				double len = MeshData::LENGTH(n);
//				if((float)len == 0){
//					flag = false;
//					break;
//				}
//				n[0] /= len;
//				n[1] /= len;
//				n[2] /= len;
//				double v1[3];
//				MeshData::VEC(v1, vertex[i], center[i1]);
//				double side1 = MeshData::DOT(n, v1);
//				double v2[3];
//				MeshData::VEC(v2, p, center[i1]);
//				double side2 = MeshData::DOT(n, v2);
//				if((float)(side1 * side2) <= 0){
//					flag = false;
//					break;
//				}
//			}
//			if(flag){
//				vertex[i][0] = (float)p[0];
//				vertex[i][1] = (float)p[1];
//				vertex[i][2] = (float)p[2];
//			}
//		}
//		else{
//		}
//	}
//
//	delete[] center;
//}

void Smoother::smoothRidgeDir(int iter, float T)
{
	int vertex_N = mesh->vertex_N;
	double (*t_max)[3] = mesh->ridge_dir;
	double (*t_min)[3] = mesh->ravine_dir;
	float (*vertex)[3] = mesh->vertex;

	double (*t_max1)[3] = new double[vertex_N][3];
	double (*t_min1)[3] = new double[vertex_N][3];

	int i;  for(i = 0; i<iter; i++){
		int j;  for(i = 0; j<vertex_N; j++){
			int *nei, nei_N;
			mesh->getConsistent1Ring(j, nei, nei_N);
			if(nei_N == 0){
				t_max1[j][0] = t_max[j][0];
				t_max1[j][1] = t_max[j][1];
				t_max1[j][2] = t_max[j][2];

				t_min1[j][0] = t_min[j][0];
				t_min1[j][1] = t_min[j][1];
				t_min1[j][2] = t_min[j][2];
				continue;
			}
			
			t_max1[j][0] = t_max1[j][1] = t_max1[j][2] = 0;
			t_min1[j][0] = t_min1[j][1] = t_min1[j][2] = 0;
			for(int k=0; k<nei_N; k++){
				int pair = nei[k];
				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);

				double dot = MeshData::DOT(t_max[j], t_max[pair]);
				double dk = acos(fabs(dot));
				double w = exp(-dk*dk/d2/(T*T));
				if(dot < 0)
					w = -w;
				t_max1[j][0] += w*t_max[pair][0];
				t_max1[j][1] += w*t_max[pair][1];
				t_max1[j][2] += w*t_max[pair][2];
				
				dot = MeshData::DOT(t_min[j], t_min[pair]);
				dk = acos(fabs(dot));
				w = exp(-dk*dk/d2/(T*T));
				if(dot < 0)
					w = -w;
				t_min1[j][0] += w*t_min[pair][0];
				t_min1[j][1] += w*t_min[pair][1];
				t_min1[j][2] += w*t_min[pair][2];
			}

			double len = MeshData::LENGTH(t_max1[j]);
			if((float)len != 0){
				t_max1[j][0] /= len;
				t_max1[j][1] /= len;
				t_max1[j][2] /= len;
			}
			else{
				t_max1[j][0] = t_max[j][0];
				t_max1[j][1] = t_max[j][1];
				t_max1[j][2] = t_max[j][2];
			}

			len = MeshData::LENGTH(t_min1[j]);
			if((float)len != 0){
				t_min1[j][0] /= len;
				t_min1[j][1] /= len;
				t_min1[j][2] /= len;
			}
			else{
				t_min1[j][0] = t_min[j][0];
				t_min1[j][1] = t_min[j][1];
				t_min1[j][2] = t_min[j][2];
			}
		}
		for(j=0; j<vertex_N; j++){
			t_max[j][0] = t_max1[j][0];
			t_max[j][1] = t_max1[j][1];
			t_max[j][2] = t_max1[j][2];
			
			t_min[j][0] = t_min1[j][0];
			t_min[j][1] = t_min1[j][1];
			t_min[j][2] = t_min1[j][2];
		}
	}
	delete[] t_max1;
	delete[] t_min1;
}

//void Smoother::smoothTmaxTmin2(int iter, float T)
//{
//	int vertex_N = mesh->vertex_N;
//	double (*t_max)[3] = mesh->t_max;
//	double (*t_min)[3] = mesh->t_min;
//	float (*vertex)[3] = mesh->vertex;
//
//	double (*t_max1)[3] = new double[vertex_N][3];
//	double (*t_min1)[3] = new double[vertex_N][3];
//
//	float (*normal)[3] = mesh->normal;
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			int *nei, nei_N;
//			mesh->getConsistent1Ring(j, nei, nei_N);
//			if(nei_N == 0){
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//				continue;
//			}
//			
//			double angle_max = 0;
//			double angle_min = 0;
//			double w_max = 0;
//			double w_min = 0;
//			for(int k=0; k<nei_N; k++){
//				int pair = nei[k];
//				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
//							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
//							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
//				
//				double cross[3];
//				MeshData::CROSS(cross, normal[pair], normal[j]);
//				double len = MeshData::LENGTH(cross);
//				double R[3][3];
//				if(len < 0.000001){
//					R[0][0] = R[1][1] = R[2][2] = 1;
//					R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
//				}
//				else{
//					double dot = MeshData::DOT(normal[pair], normal[j]);
//					if(dot > 1.0)
//						dot = 1;
//					else if(dot < -1.0)
//						dot = -1;
//					
//					GENERATE_MAT(R, acos(dot), cross);
//				}
//				double tmp1[3], tmp2[3];	
//				MAT_VEC(tmp1, R, t_max[pair]);
//				MAT_VEC(tmp2, R, t_min[pair]);
//
//				double dot = MeshData::DOT(t_max[j], tmp1);
//				if(dot > 1.0)
//					dot = 1.0;
//				else if(dot < -1.0)
//					dot = -1.0;
//				double dk = acos(fabs(dot));
//				double w = exp(-dk*dk/d2/(T*T));
//	
//				w = 1.0;
//
//				w_max += w;
//				MeshData::CROSS(cross, t_max[j], tmp1);
//				if(MeshData::DOT(cross, normal[j]) < 0)
//					dk = -dk;
//				angle_max += w*dk;
//
//				dot = MeshData::DOT(t_min[j], tmp2);
//				if(dot > 1.0)
//					dot = 1.0;
//				else if(dot < -1.0)
//					dot = -1.0;
//				dk = acos(fabs(dot));
//				w = exp(-dk*dk/d2/(T*T));
//
//				w = 1.0;
//
//				w_min += w;
//				MeshData::CROSS(cross, t_min[j], tmp2);
//				if(MeshData::DOT(cross, normal[j]) < 0)
//					dk = -dk;
//				angle_min += w*dk;
//			}
//			if(w_max != 0){
//				angle_max /= w_max;
//				double R[3][3];
//				GENERATE_MAT(R, angle_max, normal[j]);
//				MAT_VEC(t_max1[j], R, t_max[j]);
//			}
//			else{
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//			}
//			if(w_min != 0){
//				angle_min /= w_min;
//				double R[3][3];
//				GENERATE_MAT(R, angle_min, normal[j]);
//				MAT_VEC(t_min1[j], R, t_min[j]);
//			}
//			else{
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//			}
//		}
//		for(j=0; j<vertex_N; j++){
//			t_max[j][0] = t_max1[j][0];
//			t_max[j][1] = t_max1[j][1];
//			t_max[j][2] = t_max1[j][2];
//			
//			t_min[j][0] = t_min1[j][0];
//			t_min[j][1] = t_min1[j][1];
//			t_min[j][2] = t_min1[j][2];
//		}
//	}
//
//	delete[] t_max1;
//	delete[] t_min1;
//
//	for(i=0; i<vertex_N; i++){
//		if(fabs(mesh->k_max[i]) > fabs(mesh->k_min[i])){
//			double len = MeshData::LENGTH(t_max[i]);
//			if(len != 0){
//				t_max[i][0] /= len;
//				t_max[i][1] /= len;
//				t_max[i][2] /= len;
//
//				MeshData::CROSS(t_min[i], normal[i], t_max[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//		else{
//			double len = MeshData::LENGTH(t_min[i]);
//			if(len != 0){
//				t_min[i][0] /= len;
//				t_min[i][1] /= len;
//				t_min[i][2] /= len;
//
//				MeshData::CROSS(t_max[i], t_min[i], normal[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//	}
//
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, normal[i], t_max[i]);
//		if(MeshData::DOT(n, t_min[i])< 0){
//			t_min[i][0] = -t_min[i][0];
//			t_min[i][1] = -t_min[i][1];
//			t_min[i][2] = -t_min[i][2];
//		}	
//	}
//}

void Smoother::smoothRidgeDir2(int iter, float T)
{
	int vertex_N = mesh->vertex_N;
	double (*t_max)[3] = mesh->ridge_dir;
	double (*t_min)[3] = mesh->ravine_dir;
	float (*vertex)[3] = mesh->vertex;

	double (*t_max1)[3] = new double[vertex_N][3];
	double (*t_min1)[3] = new double[vertex_N][3];

	float (*normal)[3] = mesh->normal;

	int i;  for(i = 0; i<iter; i++){
		int j;  for(i = 0; j<vertex_N; j++){
			int *nei, nei_N;
			mesh->getConsistent1Ring(j, nei, nei_N);
			if(nei_N == 0){
				t_max1[j][0] = t_max[j][0];
				t_max1[j][1] = t_max[j][1];
				t_max1[j][2] = t_max[j][2];

				t_min1[j][0] = t_min[j][0];
				t_min1[j][1] = t_min[j][1];
				t_min1[j][2] = t_min[j][2];
				continue;
			}
			
			double angle_max = 0;
			double angle_min = 0;
			double w_max = 0;
			double w_min = 0;
			for(int k=0; k<nei_N; k++){
				int pair = nei[k];
				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
				
				double cross[3];
				MeshData::CROSS(cross, normal[pair], normal[j]);
				double len = MeshData::LENGTH(cross);
				double R[3][3];
				if(len < 0.000001){
					R[0][0] = R[1][1] = R[2][2] = 1;
					R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
				}
				else{
					double dot = MeshData::DOT(normal[pair], normal[j]);
					if(dot > 1.0)
						dot = 1;
					else if(dot < -1.0)
						dot = -1;
					
					GENERATE_MAT(R, acos(dot), cross);
				}
				double tmp1[3], tmp2[3];	
				MAT_VEC(tmp1, R, t_max[pair]);
				MAT_VEC(tmp2, R, t_min[pair]);

				double dot = MeshData::DOT(t_max[j], tmp1);
				if(dot > 1.0)
					dot = 1.0;
				else if(dot < -1.0)
					dot = -1.0;
				double dk = acos(fabs(dot));
				double w = exp(-dk*dk/d2/(T*T));
				w_max += w;
				MeshData::CROSS(cross, t_max[j], tmp1);
				if(MeshData::DOT(cross, normal[j]) < 0)
					dk = -dk;
				angle_max += w*dk;

				dot = MeshData::DOT(t_min[j], tmp2);
				if(dot > 1.0)
					dot = 1.0;
				else if(dot < -1.0)
					dot = -1.0;
				dk = acos(fabs(dot));
				w = exp(-dk*dk/d2/(T*T));
				w_min += w;
				MeshData::CROSS(cross, t_min[j], tmp2);
				if(MeshData::DOT(cross, normal[j]) < 0)
					dk = -dk;
				angle_min += w*dk;
			}
			if(w_max != 0){
				angle_max /= w_max;
				double R[3][3];
				GENERATE_MAT(R, angle_max, normal[j]);
				MAT_VEC(t_max1[j], R, t_max[j]);
			}
			else{
				t_max1[j][0] = t_max[j][0];
				t_max1[j][1] = t_max[j][1];
				t_max1[j][2] = t_max[j][2];
			}
			if(w_min != 0){
				angle_min /= w_min;
				double R[3][3];
				GENERATE_MAT(R, angle_min, normal[j]);
				MAT_VEC(t_min1[j], R, t_min[j]);
			}
			else{
				t_min1[j][0] = t_min[j][0];
				t_min1[j][1] = t_min[j][1];
				t_min1[j][2] = t_min[j][2];
			}
		}
		for(j=0; j<vertex_N; j++){
			t_max[j][0] = t_max1[j][0];
			t_max[j][1] = t_max1[j][1];
			t_max[j][2] = t_max1[j][2];
			
			t_min[j][0] = t_min1[j][0];
			t_min[j][1] = t_min1[j][1];
			t_min[j][2] = t_min1[j][2];
		}
	}

	delete[] t_max1;
	delete[] t_min1;

	for(i=0; i<vertex_N; i++){
		double len = MeshData::LENGTH(t_max[i]);
		if(len != 0){
			t_max[i][0] /= len;
			t_max[i][1] /= len;
			t_max[i][2] /= len;
		}
		len = MeshData::LENGTH(t_min[i]);
		if(len != 0){
			t_min[i][0] /= len;
			t_min[i][1] /= len;
			t_min[i][2] /= len;
		}
	}
}

void Smoother::smoothRidgeDir3(int iter, float T)
{
	int vertex_N = mesh->vertex_N;
	double (*t_max)[3] = mesh->ridge_dir;
	double (*t_min)[3] = mesh->ravine_dir;
	float (*vertex)[3] = mesh->vertex;

	double (*t_max1)[3] = new double[vertex_N][3];
	double (*t_min1)[3] = new double[vertex_N][3];

	float (*normal)[3] = mesh->normal;

	int i;  for(i = 0; i<iter; i++){
		int j;  for(i = 0; j<vertex_N; j++){
			int *nei, nei_N;
			mesh->getConsistent1Ring(j, nei, nei_N);
			if(nei_N == 0){
				t_max1[j][0] = t_max[j][0];
				t_max1[j][1] = t_max[j][1];
				t_max1[j][2] = t_max[j][2];

				t_min1[j][0] = t_min[j][0];
				t_min1[j][1] = t_min[j][1];
				t_min1[j][2] = t_min[j][2];
				continue;
			}
			double max_tmp[3] = {0,0,0};
			double min_tmp[3] = {0,0,0};
			for(int k=0; k<nei_N; k++){
				int pair = nei[k];
				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
				
				double cross[3];
				MeshData::CROSS(cross, normal[pair], normal[j]);
				double len = MeshData::LENGTH(cross);
				double R[3][3];
				if(len < 0.000001){
					R[0][0] = R[1][1] = R[2][2] = 1;
					R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
				}
				else{
					double dot = MeshData::DOT(normal[pair], normal[j]);
					if(dot > 1.0)
						dot = 1;
					else if(dot < -1.0)
						dot = -1;
					
					GENERATE_MAT(R, acos(dot), cross);
				}
				double tmp1[3], tmp2[3];	
				MAT_VEC(tmp1, R, t_max[pair]);
				MAT_VEC(tmp2, R, t_min[pair]);

				double dot = MeshData::DOT(t_max[j], tmp1);
				if(dot > 1.0)
					dot = 1.0;
				else if(dot < -1.0)
					dot = -1.0;
				double dk = acos(fabs(dot));
				double w = exp(-dk*dk/d2/(T*T));
				if(dot < 0)
					w = -w;
				max_tmp[0] += w*tmp1[0];
				max_tmp[1] += w*tmp1[1];
				max_tmp[2] += w*tmp1[2];

				dot = MeshData::DOT(t_min[j], tmp2);
				if(dot > 1.0)
					dot = 1.0;
				else if(dot < -1.0)
					dot = -1.0;
				dk = acos(fabs(dot));
				w = exp(-dk*dk/d2/(T*T));
				if(dot < 0)
					w = -w;
				min_tmp[0] += w*tmp2[0];
				min_tmp[1] += w*tmp2[1];
				min_tmp[2] += w*tmp2[2];
			}
			double len = MeshData::LENGTH(max_tmp);
			if((float)len != 0){
				t_max1[j][0] = max_tmp[0]/len;
				t_max1[j][1] = max_tmp[1]/len;
				t_max1[j][2] = max_tmp[2]/len;
			}
			else{
				t_max1[j][0] = t_max[j][0];
				t_max1[j][1] = t_max[j][1];
				t_max1[j][2] = t_max[j][2];
			}
			len = MeshData::LENGTH(min_tmp);
			if((float)len != 0){
				t_min1[j][0] = min_tmp[0]/len;
				t_min1[j][1] = min_tmp[1]/len;
				t_min1[j][2] = min_tmp[2]/len;
			}
			else{
				t_min1[j][0] = t_min[j][0];
				t_min1[j][1] = t_min[j][1];
				t_min1[j][2] = t_min[j][2];
			}
		}
		for(j=0; j<vertex_N; j++){
			t_max[j][0] = t_max1[j][0];
			t_max[j][1] = t_max1[j][1];
			t_max[j][2] = t_max1[j][2];
			
			t_min[j][0] = t_min1[j][0];
			t_min[j][1] = t_min1[j][1];
			t_min[j][2] = t_min1[j][2];
		}
	}

	delete[] t_max1;
	delete[] t_min1;
}

//void Smoother::smoothTmaxTmin3(int iter, float T)
//{
//	int vertex_N = mesh->vertex_N;
//	double (*t_max)[3] = mesh->t_max;
//	double (*t_min)[3] = mesh->t_min;
//	float (*vertex)[3] = mesh->vertex;
//
//	double (*t_max1)[3] = new double[vertex_N][3];
//	double (*t_min1)[3] = new double[vertex_N][3];
//
//	float (*normal)[3] = mesh->normal;
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			int *nei, nei_N;
//			mesh->getConsistent1Ring(j, nei, nei_N);
//			if(nei_N == 0){
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//				continue;
//			}
//			double max_tmp[3] = {0,0,0};
//			double min_tmp[3] = {0,0,0};
//			for(int k=0; k<nei_N; k++){
//				int pair = nei[k];
//				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
//							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
//							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
//				
//				double cross[3];
//				MeshData::CROSS(cross, normal[pair], normal[j]);
//				double len = MeshData::LENGTH(cross);
//				double R[3][3];
//				if(len < 0.000001){
//					R[0][0] = R[1][1] = R[2][2] = 1;
//					R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
//				}
//				else{
//					double dot = MeshData::DOT(normal[pair], normal[j]);
//					if(dot > 1.0)
//						dot = 1;
//					else if(dot < -1.0)
//						dot = -1;
//					
//					GENERATE_MAT(R, acos(dot), cross);
//				}
//				double tmp1[3], tmp2[3];	
//				MAT_VEC(tmp1, R, t_max[pair]);
//				MAT_VEC(tmp2, R, t_min[pair]);
//
//				double dot = MeshData::DOT(t_max[j], tmp1);
//				if(dot > 1.0)
//					dot = 1.0;
//				else if(dot < -1.0)
//					dot = -1.0;
//				double dk = acos(fabs(dot));
//				double w = exp(-dk*dk/d2/(T*T));
//
//				//w = 1.0f;
//
//				if(dot < 0)
//					w = -w;
//				max_tmp[0] += w*tmp1[0];
//				max_tmp[1] += w*tmp1[1];
//				max_tmp[2] += w*tmp1[2];
//
//				dot = MeshData::DOT(t_min[j], tmp2);
//				if(dot > 1.0)
//					dot = 1.0;
//				else if(dot < -1.0)
//					dot = -1.0;
//				dk = acos(fabs(dot));
//				w = exp(-dk*dk/d2/(T*T));
//
//				//w = 1.0f;
//
//				if(dot < 0)
//					w = -w;
//				min_tmp[0] += w*tmp2[0];
//				min_tmp[1] += w*tmp2[1];
//				min_tmp[2] += w*tmp2[2];
//			}
//			double len = MeshData::LENGTH(max_tmp);
//			if((float)len != 0){
//				t_max1[j][0] = max_tmp[0]/len;
//				t_max1[j][1] = max_tmp[1]/len;
//				t_max1[j][2] = max_tmp[2]/len;
//			}
//			else{
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//			}
//			len = MeshData::LENGTH(min_tmp);
//			if((float)len != 0){
//				t_min1[j][0] = min_tmp[0]/len;
//				t_min1[j][1] = min_tmp[1]/len;
//				t_min1[j][2] = min_tmp[2]/len;
//			}
//			else{
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//			}
//		}
//		for(j=0; j<vertex_N; j++){
//			t_max[j][0] = t_max1[j][0];
//			t_max[j][1] = t_max1[j][1];
//			t_max[j][2] = t_max1[j][2];
//			
//			t_min[j][0] = t_min1[j][0];
//			t_min[j][1] = t_min1[j][1];
//			t_min[j][2] = t_min1[j][2];
//		}
//	}
//
//	delete[] t_max1;
//	delete[] t_min1;
//
//	for(i=0; i<vertex_N; i++){
//		if(fabs(mesh->k_max[i]) > fabs(mesh->k_min[i])){
//			double len = MeshData::LENGTH(t_max[i]);
//			if(len != 0){
//				t_max[i][0] /= len;
//				t_max[i][1] /= len;
//				t_max[i][2] /= len;
//
//				MeshData::CROSS(t_min[i], normal[i], t_max[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//		else{
//			double len = MeshData::LENGTH(t_min[i]);
//			if(len != 0){
//				t_min[i][0] /= len;
//				t_min[i][1] /= len;
//				t_min[i][2] /= len;
//
//				MeshData::CROSS(t_max[i], t_min[i], normal[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//	}
//
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, normal[i], t_max[i]);
//		if(MeshData::DOT(n, t_min[i])< 0){
//			t_min[i][0] = -t_min[i][0];
//			t_min[i][1] = -t_min[i][1];
//			t_min[i][2] = -t_min[i][2];
//		}	
//	}
//}

void Smoother::smoothKmaxKmin2(int iter, float T)
{
	int vertex_N = mesh->vertex_N;
	double *k_max = mesh->k_max;
	double *k_min = mesh->k_min;
	float (*vertex)[3] = mesh->vertex;

	double *k_max1 = new double[vertex_N];
	double *k_min1 = new double[vertex_N];

	int i;  for(i = 0; i<iter; i++){
		float dt = T;
		//if(i%2 != 0)
			//dt = (float)(dt/(dt*T - 1.0));

		int j;  for(i = 0; j<vertex_N; j++){
			int *nei, nei_N;

			k_max1[j] = 0;//k_max[j];
			k_min1[j] = 0;//k_min[j];

			mesh->getConsistent1Ring(j, nei, nei_N);
			if(nei_N == 0){
				continue;
			}
			
			double total_max = 0;//1;
			double total_min = 0;//1;
		
			for(int k=0; k<nei_N; k++){
				int pair = nei[k];
				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);

				//double dk = k_max[j] - k_max[pair];
				//double w = exp(-dk*dk/d2/(T*T));
				//float w = exp(-d2/(T*T));
				float w = 1.0;
				total_max += w;
				k_max1[j] += w*k_max[pair];

				//dk = k_min[j] - k_min[pair];
				//w = exp(-dk*dk/d2/(T*T));
				//w = exp(-d2/(T*T));
				w = 1.0;
				total_min += w;
				k_min1[j] += w*k_min[pair];
			}
			if(total_max != 0)
				k_max1[j] = (1.0-dt)*k_max[j] + dt*k_max1[j]/total_max;
			else
				k_max1[j] = k_max[j];

			if(total_min != 0)
				k_min1[j] = (1.0-dt)*k_min[j] + dt*k_min1[j]/total_min;
				//k_min1[j] /= total_min;
			else
				k_min1[j] = k_min[j];
		}
		for(j=0; j<vertex_N; j++){
			k_max[j] = k_max1[j];
			k_min[j] = k_min1[j];
		}
	}

	delete[] k_max1;
	delete[] k_min1;
}

//void Smoother::smoothTmaxTmin4(int iter, float T)
//{
//	int vertex_N = mesh->vertex_N;
//	double (*t_max)[3] = mesh->t_max;
//	double (*t_min)[3] = mesh->t_min;
//	float (*vertex)[3] = mesh->vertex;
//
//	double (*t_max1)[3] = new double[vertex_N][3];
//	double (*t_min1)[3] = new double[vertex_N][3];
//
//	float (*normal)[3] = mesh->normal;
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			int *nei, nei_N;
//			mesh->getConsistent1Ring(j, nei, nei_N);
//			if(nei_N == 0){
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//				continue;
//			}
//			double max_tmp[3];
//			max_tmp[0] = t_max[j][0];
//			max_tmp[1] = t_max[j][1];
//			max_tmp[2] = t_max[j][2];
//			double min_tmp[3];
//			min_tmp[0] = t_min[j][0];
//			min_tmp[1] = t_min[j][1];
//			min_tmp[2] = t_min[j][2];
//
//			for(int k=0; k<nei_N; k++){
//				int pair = nei[k];
//				double d2 = (vertex[j][0]-vertex[pair][0])*(vertex[j][0]-vertex[pair][0])
//							+ (vertex[j][1]-vertex[pair][1])*(vertex[j][1]-vertex[pair][1])
//							+ (vertex[j][2]-vertex[pair][2])*(vertex[j][2]-vertex[pair][2]);
//				
//				double cross[3];
//				MeshData::CROSS(cross, normal[pair], normal[j]);
//				double len = MeshData::LENGTH(cross);
//				double R[3][3];
//				if(len < 0.000001){
//					R[0][0] = R[1][1] = R[2][2] = 1;
//					R[0][1] = R[0][2] = R[1][0] = R[1][2] = R[2][0] = R[2][1] = 0;
//				}
//				else{
//					double dot = MeshData::DOT(normal[pair], normal[j]);
//					if(dot > 1.0)
//						dot = 1;
//					else if(dot < -1.0)
//						dot = -1;
//					
//					GENERATE_MAT(R, acos(dot), cross);
//				}
//				double tmp1[3], tmp2[3];	
//				MAT_VEC(tmp1, R, t_max[pair]);
//				MAT_VEC(tmp2, R, t_min[pair]);
//
//				double dot = MeshData::DOT(t_max[j], tmp1);
//				if(dot > 1.0)
//					dot = 1.0;
//				else if(dot < -1.0)
//					dot = -1.0;
//				double dk = acos(fabs(dot));
//				double w = exp(-dk*dk/d2/(T*T));
//
//				//w = 1.0f;
//
//				if(dot < 0)
//					w = -w;
//				max_tmp[0] += w*tmp1[0];
//				max_tmp[1] += w*tmp1[1];
//				max_tmp[2] += w*tmp1[2];
//
//				dot = MeshData::DOT(t_min[j], tmp2);
//				if(dot > 1.0)
//					dot = 1.0;
//				else if(dot < -1.0)
//					dot = -1.0;
//				dk = acos(fabs(dot));
//				w = exp(-dk*dk/d2/(T*T));
//
//				//w = 1.0f;
//
//				if(dot < 0)
//					w = -w;
//				min_tmp[0] += w*tmp2[0];
//				min_tmp[1] += w*tmp2[1];
//				min_tmp[2] += w*tmp2[2];
//			}
//			double len = MeshData::LENGTH(max_tmp);
//			if((float)len != 0){
//				t_max1[j][0] = max_tmp[0]/len;
//				t_max1[j][1] = max_tmp[1]/len;
//				t_max1[j][2] = max_tmp[2]/len;
//			}
//			else{
//				t_max1[j][0] = t_max[j][0];
//				t_max1[j][1] = t_max[j][1];
//				t_max1[j][2] = t_max[j][2];
//			}
//			len = MeshData::LENGTH(min_tmp);
//			if((float)len != 0){
//				t_min1[j][0] = min_tmp[0]/len;
//				t_min1[j][1] = min_tmp[1]/len;
//				t_min1[j][2] = min_tmp[2]/len;
//			}
//			else{
//				t_min1[j][0] = t_min[j][0];
//				t_min1[j][1] = t_min[j][1];
//				t_min1[j][2] = t_min[j][2];
//			}
//		}
//		for(j=0; j<vertex_N; j++){
//			t_max[j][0] = t_max1[j][0];
//			t_max[j][1] = t_max1[j][1];
//			t_max[j][2] = t_max1[j][2];
//			
//			t_min[j][0] = t_min1[j][0];
//			t_min[j][1] = t_min1[j][1];
//			t_min[j][2] = t_min1[j][2];
//		}
//	}
//
//	delete[] t_max1;
//	delete[] t_min1;
//
//	for(i=0; i<vertex_N; i++){
//		if(fabs(mesh->k_max[i]) > fabs(mesh->k_min[i])){
//			double len = MeshData::LENGTH(t_max[i]);
//			if(len != 0){
//				t_max[i][0] /= len;
//				t_max[i][1] /= len;
//				t_max[i][2] /= len;
//
//				MeshData::CROSS(t_min[i], normal[i], t_max[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//		else{
//			double len = MeshData::LENGTH(t_min[i]);
//			if(len != 0){
//				t_min[i][0] /= len;
//				t_min[i][1] /= len;
//				t_min[i][2] /= len;
//
//				MeshData::CROSS(t_max[i], t_min[i], normal[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//	}
//
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, normal[i], t_max[i]);
//		if(MeshData::DOT(n, t_min[i])< 0){
//			t_min[i][0] = -t_min[i][0];
//			t_min[i][1] = -t_min[i][1];
//			t_min[i][2] = -t_min[i][2];
//		}	
//	}
//}

//void Smoother::smoothNormalL()
//{
//	float (*normal)[3] = mesh->normal_f;
//	int (*face)[3] = mesh->face;
//	int face_N = mesh->face_N;
//	float (*tmp_normal)[3] = new float[face_N][3];
//	int (*ad_face)[3] = mesh->face_link_E;
//
//	int i;  for(i = 0; i<face_N; i++){
//		double sum_n[3];
//		sum_n[0] = normal[i][0];
//		sum_n[1] = normal[i][1];
//		sum_n[2] = normal[i][2];
//		int j;  for(i = 0; j<3; j++){
//			int pair = ad_face[i][j];
//			if(pair < 0)
//				continue;
//			/*
//double c1[3], c2[3];
//mesh->faceCenter(c1, i);
//mesh->faceCenter(c2, pair);
//double d = MeshData::DIST(c1, c2);
//double a = 0;
//if(d != 0)
//	a = 1.0/d;
//	*/
//			double a = 1;
//			sum_n[0] += a*normal[pair][0];
//			sum_n[1] += a*normal[pair][1];
//			sum_n[2] += a*normal[pair][2];
//		}
//		double len = MeshData::LENGTH(sum_n);
//		if(len != 0){
//			tmp_normal[i][0] = (float)(sum_n[0]/len);
//			tmp_normal[i][1] = (float)(sum_n[1]/len);
//			tmp_normal[i][2] = (float)(sum_n[2]/len);
//		}
//		else{
//			tmp_normal[i][0] = normal[i][0];
//			tmp_normal[i][1] = normal[i][1];
//			tmp_normal[i][2] = normal[i][2];
//		}
//	}
//	delete[] normal;
//	mesh->normal_f = tmp_normal;
//}

//void Smoother::smoothTmaxTminSIG1(int iter, float T)
//T is not used.
//{
//	int vertex_N = mesh->vertex_N;
//	double (*t_max)[3] = mesh->t_max;
//	double (*t_min)[3] = mesh->t_min;
//	float (*vertex)[3] = mesh->vertex;
//	int *degree = mesh->degree_v;
//	int **link = mesh->vertex_link_v;
//	float (*normal)[3] = mesh->normal;
//
//	double *theta_max = new double[vertex_N];
//	double *theta_min = new double[vertex_N];
//
//	double **phi = new double*[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++){
//		int deg = degree[i];
//		if(deg == 0)
//			continue;
//		int *l = link[i]; 
//		phi[i] = new double[deg];
//		double* p = phi[i];
//		p[0] = 0;
//		double t0[3];
//		mesh->tangent(t0, i, l[0]); 
//		for(int j=1; j<deg; j++){
//			double tj[3];
//			mesh->tangent(tj, i, l[j]);
//			double dot = MeshData::DOT(t0, tj);
//			if(dot > 1)
//				dot = 1;
//			else if(dot < -1)
//				dot = -1;
//			p[j] = acos(dot);
//			double c[3];
//			MeshData::CROSS(c, t0, tj);
//			if(MeshData::DOT(c, normal[i]) < 0)
//				p[j] = -p[j];
//		}
//		double dot_max = MeshData::DOT(t0, t_max[i]);
//		if(dot_max > 1)
//			dot_max = 1;
//		else if(dot_max < -1)
//			dot_max = -1;
//		theta_max[i] = acos(dot_max);
//		double c[3];
//		MeshData::CROSS(c, t0, t_max[i]);
//		if(MeshData::DOT(c, normal[i]) < 0)
//			theta_max[i] = -theta_max[i];
//
//		double dot_min = MeshData::DOT(t0, t_min[i]);
//		if(dot_min > 1)
//			dot_min = 1;
//		else if(dot_min < -1)
//			dot_min = -1;
//		theta_min[i] = acos(dot_min);
//		MeshData::CROSS(c, t0, t_min[i]);
//		if(MeshData::DOT(c, normal[i]) < 0)
//			theta_min[i] = -theta_min[i];
//	}
//
//	double *theta_max1 = new double[vertex_N];
//	double *theta_min1 = new double[vertex_N];
//
//	for(i=0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			int deg = degree[j];
//			if(deg == 0){
//				theta_max1[j] = theta_max[j];
//				theta_min1[j] = theta_min[j];
//			}
//			double* p = phi[j];
//			int* l = link[j];
//			double d_theta_max = 0;
//			double d_theta_min = 0;
//			for(int k=0; k<deg; k++){
//				int pair = l[k];
//				int search = -1;
//				for(int m=0; m<degree[pair]; m++)
//					if(link[pair][m] == j){
//						search = m;
//						break;
//					}
//				if(search < 0)
//					continue;
//				d_theta_max += 2.0*sin(2.0*((theta_max[j]-p[k]) 
//					- (theta_max[pair]-phi[pair][search])));
//				d_theta_min += 2.0*sin(2.0*((theta_min[j]-p[k]) 
//					- (theta_min[pair]-phi[pair][search])));	
//			}
//			theta_max1[j] = theta_max[j] - 0.1*d_theta_max/(double)deg;
//			theta_min1[j] = theta_min[j] - 0.1*d_theta_min/(double)deg;
//		}
//	}
//	
//	// compute t_{max,min} from theta_{max,min}
//	for(i=0; i<vertex_N; i++){
//		double R[3][3];
//		double axis[3];
//		axis[0] = (double)normal[i][0];
//		axis[1] = (double)normal[i][1];
//		axis[2] = (double)normal[i][2];
//		double t0[3];
//		mesh->tangent(t0, i, link[i][0]);
//
//		GENERATE_MAT(R, theta_max[i], axis);
//		MAT_VEC(t_max[i], R, t0);
//
//		GENERATE_MAT(R, theta_min[i], axis);
//		MAT_VEC(t_min[i], R, t0);
//	}
//
// 	delete[] theta_max;
//	delete[] theta_min;
//	for(i=0; i<vertex_N; i++)
//		if(degree[i] != 0)
//			delete[] phi[i];
//	delete[] phi;
//
//	for(i=0; i<vertex_N; i++){
//		if(fabs(mesh->k_max[i]) > fabs(mesh->k_min[i])){
//			double len = MeshData::LENGTH(t_max[i]);
//			if(len != 0){
//				t_max[i][0] /= len;
//				t_max[i][1] /= len;
//				t_max[i][2] /= len;
//
//				MeshData::CROSS(t_min[i], normal[i], t_max[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//		else{
//			double len = MeshData::LENGTH(t_min[i]);
//			if(len != 0){
//				t_min[i][0] /= len;
//				t_min[i][1] /= len;
//				t_min[i][2] /= len;
//
//				MeshData::CROSS(t_max[i], t_min[i], normal[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//	}
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, normal[i], t_max[i]);
//		if(MeshData::DOT(n, t_min[i])< 0){
//			t_min[i][0] = -t_min[i][0];
//			t_min[i][1] = -t_min[i][1];
//			t_min[i][2] = -t_min[i][2];
//		}	
//	}
//}

//void Smoother::smoothTmaxTminSIG2()
//{
//	int vertex_N = mesh->vertex_N;
//	double (*t_max)[3] = mesh->t_max;
//	double (*t_min)[3] = mesh->t_min;
//	float (*vertex)[3] = mesh->vertex;
//	int *degree = mesh->degree_v;
//	int **link = mesh->vertex_link_v;
//	float (*normal)[3] = mesh->normal;
//
//	float *theta_max = new float[vertex_N];
//	float *theta_min = new float[vertex_N];
//
//	double **phi = new double*[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++){
//		int deg = degree[i];
//		if(deg == 0)
//			continue;
//		int *l = link[i]; 
//		phi[i] = new double[deg];
//		double* p = phi[i];
//		p[0] = 0;
//		double t0[3];
//		mesh->tangent(t0, i, l[0]); 
//		for(int j=1; j<deg; j++){
//			double tj[3];
//			mesh->tangent(tj, i, l[j]);
//			double dot = MeshData::DOT(t0, tj);
//			if(dot > 1)
//				dot = 1;
//			else if(dot < -1)
//				dot = -1;
//			p[j] = acos(dot);
//			double c[3];
//			MeshData::CROSS(c, t0, tj);
//			if(MeshData::DOT(c, normal[i]) < 0)
//				p[j] = -p[j];
//		}
//		double dot_max = MeshData::DOT(t0, t_max[i]);
//		if(dot_max > 1)
//			dot_max = 1;
//		else if(dot_max < -1)
//			dot_max = -1;
//		theta_max[i] = (float)acos(dot_max);
//		double c[3];
//		MeshData::CROSS(c, t0, t_max[i]);
//		if(MeshData::DOT(c, normal[i]) < 0)
//			theta_max[i] = -theta_max[i];
//
//		double dot_min = MeshData::DOT(t0, t_min[i]);
//		if(dot_min > 1)
//			dot_min = 1;
//		else if(dot_min < -1)
//			dot_min = -1;
//		theta_min[i] = (float)acos(dot_min);
//		MeshData::CROSS(c, t0, t_min[i]);
//		if(MeshData::DOT(c, normal[i]) < 0)
//			theta_min[i] = -theta_min[i];
//	}
//
//	//conjugate gradient algorithm
//	FieldEnergy e; // = new FieldEnergy;
//	e.initData(vertex_N, phi, link, degree);
//	EnergyMinimizer mini; // = new EnergyMinimizer;
//	mini.energy = &e;
//	int iter;
//	float fret;
//	for(i=0; i<5; i++){
//		mini.frprmn(theta_max, vertex_N, 0.0001f, &iter, &fret);
//		mini.frprmn(theta_min, vertex_N, 0.0001f, &iter, &fret);
//	}
//
//	// compute t_{max,min} from theta_{max,min}
//	for(i=0; i<vertex_N; i++){
//		double R[3][3];
//		double axis[3];
//		axis[0] = (double)normal[i][0];
//		axis[1] = (double)normal[i][1];
//		axis[2] = (double)normal[i][2];
//		double t0[3];
//		mesh->tangent(t0, i, link[i][0]);
//
//		GENERATE_MAT(R, theta_max[i], axis);
//		MAT_VEC(t_max[i], R, t0);
//
//		GENERATE_MAT(R, theta_min[i], axis);
//		MAT_VEC(t_min[i], R, t0);
//	}
//
// 	delete[] theta_max;
//	delete[] theta_min;
//	for(i=0; i<vertex_N; i++)
//		if(degree[i] != 0)
//			delete[] phi[i];
//	delete[] phi;
//
//	for(i=0; i<vertex_N; i++){
//		if(fabs(mesh->k_max[i]) > fabs(mesh->k_min[i])){
//			double len = MeshData::LENGTH(t_max[i]);
//			if(len != 0){
//				t_max[i][0] /= len;
//				t_max[i][1] /= len;
//				t_max[i][2] /= len;
//
//				MeshData::CROSS(t_min[i], normal[i], t_max[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//		else{
//			double len = MeshData::LENGTH(t_min[i]);
//			if(len != 0){
//				t_min[i][0] /= len;
//				t_min[i][1] /= len;
//				t_min[i][2] /= len;
//
//				MeshData::CROSS(t_max[i], t_min[i], normal[i]);
//			}
//			else{
//				t_max[i][0] = t_max[i][1] = t_max[i][2] = 0;
//				t_min[i][0] = t_min[i][1] = t_min[i][2] = 0;
//			}
//		}
//	}
//
//	for(i=0; i<vertex_N; i++){
//		double n[3];
//		MeshData::CROSS(n, normal[i], t_max[i]);
//		if(MeshData::DOT(n, t_min[i])< 0){
//			t_min[i][0] = -t_min[i][0];
//			t_min[i][1] = -t_min[i][1];
//			t_min[i][2] = -t_min[i][2];
//		}	
//	}
//}

//void Smoother::moveNormalD(float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	int (*face)[3] = mesh->face;
//	float (*normal_f)[3] = mesh->normal_f;
//	double (*dp)[3] = new double[vertex_N][3];
//	int **link_f = mesh->vertex_link_f;
//	int *degree = mesh->degree_v;
//
//	int i;  for(i = 0; i<vertex_N; i++){
//		if(mesh->isBound[i]){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//
//		int *nei = link_f[i];
//		int nei_N = degree[i];
//
//		if(nei_N == 0){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//
//		double total_w = 0;
//
//		double df[3];
//		df[0] = df[1] = df[2] = 0;
//
//		for(int j = 0; j<nei_N; j++){
//			int f = nei[j];
//
//			int i0 = i;
//			int i1, i2;
//			if(i0 == face[f][0]){
//				i1 = face[f][1];
//				i2 = face[f][2];
//			}
//			else if(i0 == face[f][1]){
//				i1 = face[f][2];
//				i2 = face[f][0];
//			}
//			else{
//				i1 = face[f][0];
//				i2 = face[f][1];
//			}
//
//			float v1[3], v2[3];
//			MeshData::VEC(v1, vertex[i0], vertex[i1]);
//			MeshData::VEC(v2, vertex[i0], vertex[i2]);
//			double c[3];
//			MeshData::CROSS(c, v1, v2);
//
//			double cx[3], cy[3], cz[3];
//
//			cx[0] = 0;
//			cx[1] = v1[2] - v2[2];
//			cx[2] = v2[1] - v1[1];
//			
//			cy[0] = v2[2] - v1[2];
//			cy[1] = 0;
//			cy[2] = v1[0] - v2[0];
//
//			cz[0] = v1[1] - v2[1];
//			cz[1] = v2[0] - v1[0];
//			cz[2] = 0;
//	
//			double m[3];
//			m[0] = normal_f[f][0];
//			m[1] = normal_f[f][1];
//			m[2] = normal_f[f][2];
//
//			double c1 = MeshData::LENGTH(c);
//			double nm[3];
//			nm[0] = c[0]/c1 - m[0];
//			nm[1] = c[1]/c1 - m[1];
//			nm[2] = c[2]/c1 - m[2];
//
//			double v[3];
//			v[0] = -MeshData::DOT(cx, nm);
//			v[1] = -MeshData::DOT(cy, nm);
//			v[2] = -MeshData::DOT(cz, nm);
//			
//			/*
//			double nx[3], ny[3], nz[3], dot;
//			double c2 = c[0]*c[0] + c[1]*c[1] + c[2]*c[2];
//			
//			dot = MeshData::DOT(c, cx);
//			nx[0] = c2*cx[0] - dot*c[0];
//			nx[1] = c2*cx[1] - dot*c[1];
//			nx[2] = c2*cx[2] - dot*c[2];
//
//			dot = MeshData::DOT(c, cy);
//			ny[0] = c2*cy[0] - dot*c[0];
//			ny[1] = c2*cy[1] - dot*c[1];
//			ny[2] = c2*cy[2] - dot*c[2];
//
//			dot = MeshData::DOT(c, cz);
//			nz[0] = c2*cz[0] - dot*c[0];
//			nz[1] = c2*cz[1] - dot*c[1];
//			nz[2] = c2*cz[2] - dot*c[2];
//
//			double m[3];
//			m[0] = normal_f[f][0];
//			m[1] = normal_f[f][1];
//			m[2] = normal_f[f][2];
//
//			double v[3];
//			double c3 = pow(c2, 1.5);
//			v[0] = MeshData::DOT(nx, m)/c3;
//			v[1] = MeshData::DOT(ny, m)/c3;
//			v[2] = MeshData::DOT(nz, m)/c3;
//			*/
//			double w = 1; //sqrt(c2);
//			total_w += w;
//	
//			df[0] += w*v[0];
//			df[1] += w*v[1];
//			df[2] += w*v[2];
//		}
//		if(total_w != 0){
//			dp[i][0] = df[0]/total_w;
//			dp[i][1] = df[1]/total_w;
//			dp[i][2] = df[2]/total_w;
//		}
//		else
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//	}
//
//	for(i=0; i<vertex_N; i++){
//		vertex[i][0] -= dt*(float)dp[i][0];
//		vertex[i][1] -= dt*(float)dp[i][1];
//		vertex[i][2] -= dt*(float)dp[i][2];
//	}
//	delete[] dp;
//}

//void Smoother::BilaplacianWithRidge(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	BOOL *isFeature = new BOOL[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++)
//		isFeature[i] = (mesh->isRidge[i] || mesh->isRavine[i]);
//
//	//L is Lplacian vectors
//	double (*L)[3] = new double[vertex_N][3];
//	double (*LL)[3] = new double[vertex_N][3];
//	int **link = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//	BOOL *isBound = mesh->isBound;
//
//	for(i=0; i<iter; i++){
//		//mesh->computeNormal();
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(isBound[j]){
//				L[j][0] = L[j][1] = L[j][2] = 0;
//				continue;
//			}
//
//			if(isFeature[j]){
//				L[j][0] = L[j][1] = L[j][2] = 0;
//				int *l = link[j];
//				int deg = degree[j];
//				int count = 0;
//				for(int k=0; k<deg; k++){
//					if(isFeature[l[k]]){
//						count++;
//						L[j][0] += vertex[l[k]][0];
//						L[j][1] += vertex[l[k]][1];
//						L[j][2] += vertex[l[k]][2];
//					}
//				}
//				if(count != 0){
//					L[j][0] = L[j][0]/count - vertex[j][0];
//					L[j][1] = L[j][1]/count - vertex[j][1];
//					L[j][2] = L[j][2]/count - vertex[j][2];
//				}
//				else
//					mesh->laplacian(j, L[j]);
//			}
//			else
//				mesh->laplacian(j, L[j]);
//		}
//
//		for(j=0; j<vertex_N; j++){
//			LL[j][0] = LL[j][1] = LL[j][2] = 0;
//			if(isBound[j])
//				continue;
//
//			int *l = link[j];
//			int deg = degree[j];
//			for(int k=0; k<deg; k++){
//				LL[j][0] += L[l[k]][0];
//				LL[j][1] += L[l[k]][1];
//				LL[j][2] += L[l[k]][2];
//			}
//			LL[j][0] = LL[j][0]/deg - L[j][0];
//			LL[j][1] = LL[j][1]/deg - L[j][1];
//			LL[j][2] = LL[j][2]/deg - L[j][2];
//		
//			/*
//			if(cnt == 1){
//				mesh->laplacian(j, L[j]);
//				double inner = MeshData::DOT(mesh->normal[j], L[j]);
//				L[i][0] = inner*mesh->normal[j][0];
//				L[i][1] = inner*mesh->normal[j][1];
//				L[i][2] = inner*mesh->normal[j][2];
//			}*/
//
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			float *p = mesh->vertex[j];
//			//replace each vertex
//			p[0] -= (float)(dt*LL[j][0]);
//			p[1] -= (float)(dt*LL[j][1]);
//			p[2] -= (float)(dt*LL[j][2]);
//		}
//	}
//	delete[] L;
//	delete[] LL;
//	delete[] isFeature;
//}

//void Smoother::BilaplacianWithRidge2(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	BOOL *isFeature = new BOOL[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++)
//		isFeature[i] = (mesh->isRidge[i] || mesh->isRavine[i]);
//
//	//L is Lplacian vectors
//	double (*L)[3] = new double[vertex_N][3];
//	double (*LL)[3] = new double[vertex_N][3];
//	int **link = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//	BOOL *isBound = mesh->isBound;
//
//	for(i=0; i<iter; i++){
//		mesh->computeNormal();
//		int j;  for(i = 0; j<vertex_N; j++)
//			mesh->laplacian(j, L[j]);
//
//		for(j=0; j<vertex_N; j++){
//			LL[j][0] = LL[j][1] = LL[j][2] = 0;
//			if(isBound[j])
//				continue;
//
//			int *l = link[j];
//			int deg = degree[j];
//			int count = 0;
//			BOOL f = isFeature[j];
//			for(int k=0; k<deg; k++){
//				int ad = l[k];
//				if(f && isFeature[ad] || !f){
//					count++;
//					LL[j][0] += L[ad][0];
//					LL[j][1] += L[ad][1];
//					LL[j][2] += L[ad][2];
//				}
//			}
//			if(count != 0){
//				LL[j][0] = LL[j][0]/count - L[j][0];
//				LL[j][1] = LL[j][1]/count - L[j][1];
//				LL[j][2] = LL[j][2]/count - L[j][2];
//			}
//			/*
//			if(cnt == 1){
//				mesh->laplacian(j, L[j]);
//				double inner = MeshData::DOT(mesh->normal[j], L[j]);
//				L[i][0] = inner*mesh->normal[j][0];
//				L[i][1] = inner*mesh->normal[j][1];
//				L[i][2] = inner*mesh->normal[j][2];
//			}*/
//
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			float *p = vertex[j];
//			//replace each vertex
//			p[0] -= (float)(dt*LL[j][0]);
//			p[1] -= (float)(dt*LL[j][1]);
//			p[2] -= (float)(dt*LL[j][2]);
//		}
//	}
//	delete[] L;
//	delete[] LL;
//	delete[] isFeature;
//}

//void Smoother::moveNormalTaubin(float C)
//{
//	int vertex_N = mesh->vertex_N;
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	float (*vertex)[3] = mesh->vertex;
//	int (*face)[3] = mesh->face;
//	double (*dp)[3] = new double[vertex_N][3];
//	double *total_w = new double[vertex_N];
//
//	int i;  for(i = 0; i<vertex_N; i++){
//		dp[i][0] = dp[i][1] = dp[i][2] = 0;
//		total_w[i] = 0;
//	}
//
//	for(i=0; i<face_N; i++){
//		int i0 = face[i][0];
//		int i1 = face[i][1];
//		int i2 = face[i][2];
//		float *n = normal[i];
//
//		float v1[3], v2[3];
//		MeshData::VEC(v1, vertex[i0], vertex[i1]);
//		MeshData::VEC(v2, vertex[i0], vertex[i2]);
//		double w1 = fabs(MeshData::DOT(n, v1));
//		double w2 = fabs(MeshData::DOT(n, v2));
//		dp[i0][0] += w1*v1[0] + w2*v2[0];
//		dp[i0][1] += w1*v1[1] + w2*v2[1];
//		dp[i0][2] += w1*v1[2] + w2*v2[2];
//		total_w[i0] += w1 + w2;
//
//		MeshData::VEC(v1, vertex[i1], vertex[i2]);
//		MeshData::VEC(v2, vertex[i1], vertex[i0]);
//		w1 = fabs(MeshData::DOT(n, v1));
//		w2 = fabs(MeshData::DOT(n, v2));
//		dp[i1][0] += w1*v1[0] + w2*v2[0];
//		dp[i1][1] += w1*v1[1] + w2*v2[1];
//		dp[i1][2] += w1*v1[2] + w2*v2[2];
//		total_w[i1] += w1 + w2;
//
//		MeshData::VEC(v1, vertex[i2], vertex[i0]);
//		MeshData::VEC(v2, vertex[i2], vertex[i1]);
//		w1 = fabs(MeshData::DOT(n, v1));
//		w2 = fabs(MeshData::DOT(n, v2));
//		dp[i2][0] += w1*v1[0] + w2*v2[0];
//		dp[i2][1] += w1*v1[1] + w2*v2[1];
//		dp[i2][2] += w1*v1[2] + w2*v2[2];
//		total_w[i2] += w1 + w2;
//	}
//
//	for(i=0; i<vertex_N; i++){
//		//if(fabs(total_w[i]) < 0.001)
//			//continue;
//
//		vertex[i][0] += (float)(C*dp[i][0]/total_w[i]);
//		vertex[i][1] += (float)(C*dp[i][1]/total_w[i]);
//		vertex[i][2] += (float)(C*dp[i][2]/total_w[i]);
//	}
//
//	delete[] dp;
//	delete[] total_w;
//}

//void Smoother::LaplacianLocalControl(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	double (*M)[3] = new double[vertex_N][3];
//	double (*L)[3] = new double[vertex_N][3];
//	int **link = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//
//	int i;  for(i = 0; i<vertex_N; i++)
//		mesh->laplacian(i, M[i]);
//		//mesh->meanCurvatureNormal(i, M[i]);
//
//	for(int k=0; k<10; k++){
//		for(i=0; i<vertex_N; i++){
//			L[i][0] = L[i][1] = L[i][2] = 0;
//			int deg = degree[i];
//			int *l = link[i];
//			int j;  for(i = 0; j<deg; j++){
//				int pair = l[j];
//				L[i][0] += M[pair][0];
//				L[i][1] += M[pair][1];
//				L[i][2] += M[pair][2];
//			}
//			L[i][0] = 0.9*M[i][0] + 0.1*L[i][0]/deg;
//			L[i][1] = 0.9*M[i][1] + 0.1*L[i][1]/deg;
//			L[i][2] = 0.9*M[i][2] + 0.1*L[i][2]/deg;
//		}
//		for(i=0; i<vertex_N; i++){
//			M[i][0] = L[i][0];
//			M[i][1] = L[i][1];
//			M[i][2] = L[i][2];
//		}
//	}
//
//	for(k=0; k<iter; k++){
//		mesh->computeNormal();
//		float (*normal)[3] = mesh->normal;
//		for(i=0; i<vertex_N; i++){
//			mesh->laplacian(i, L[i]);
//			/*
//			mesh->meanCurvatureNormal(i, L[i]);
//
//			double U[3];
//			float *n = normal[i];
//			mesh->laplacian(i, U);
//			double dot = MeshData::DOT(U, n);
//			U[0] = U[0] - dot*n[0];
//			U[1] = U[1] - dot*n[1];
//			U[2] = U[2] - dot*n[2];
//
//			L[i][0] -= 0.1*U[0]/dt;
//			L[i][1] -= 0.1*U[1]/dt;
//			L[i][2] -= 0.1*U[2]/dt;
//			*/
//		}
//
//		for(i=0; i<vertex_N; i++){
//			//vertex[i][0] -= (float)(dt*L[i][0]);
//			//vertex[i][1] -= (float)(dt*L[i][1]);
//			//vertex[i][2] -= (float)(dt*L[i][2]);
//			vertex[i][0] += (float)(dt*(L[i][0] - M[i][0]));
//			vertex[i][1] += (float)(dt*(L[i][1] - M[i][1]));
//			vertex[i][2] += (float)(dt*(L[i][2] - M[i][2]));
//		}
//	}
//}

void Smoother::TaubinMethod(int iter, float dt, float low_pass, int type)
{
	int vertex_N = mesh->vertex_N;
	int *degree = mesh->degree_v;
	int **link = mesh->vertex_link_v;
	float (*vertex)[3] = mesh->vertex;
	double (*L)[3] = new double[vertex_N][3];

	float dt1 = dt;
	float dt2 = (float)(dt/(dt*low_pass - 1.0));

	int i;  for(i = 0; i<2*iter; i++){
		if(i%2 == 0)
			dt = dt1;
		else
			dt = dt2;

		if(type == 0){
			int j;  for(i = 0; j<vertex_N; j++){
				//Boundary points are fixed
				if(mesh->isBound[j]){
					L[j][0] = L[j][1] = L[j][2] = 0;
				}	
				else{
					mesh->laplacian(j, L[j]);
				}
			}
		}
		else if(type == 1){
			int j;  for(i = 0; j<vertex_N; j++){
				L[j][0] = L[j][1] = L[j][2] = 0;

				//Boundary points are fixed
				if(mesh->isBound[j])
					continue;
	
				double total_w = 0;
				int deg = degree[j];
				int *l = link[j];
				float *p = vertex[j];
				for(int k=0; k<deg; k++){
					float *p1 = vertex[l[k]];
					double w = MeshData::DIST(p, p1);
					if((float)w == 0)
						continue;
					w = 1.0/w;
					total_w += w;
					float v[3];
					MeshData::VEC(v, p, p1);
					L[j][0] += w*v[0];
					L[j][1] += w*v[1];
					L[j][2] += w*v[2];
				}
				if(total_w != 0){
					L[j][0] /= total_w;
					L[j][1] /= total_w;
					L[j][2] /= total_w;
				}
				else
					L[j][0] = L[j][1] = L[j][2] = 0;
			}
		}
		else{
			int j;  for(i = 0; j<vertex_N; j++){
				L[j][0] = L[j][1] = L[j][2] = 0;

				//Boundary points are fixed
				if(mesh->isBound[j])
					continue;
	
				double total_w = 0;
				int deg = degree[j];
				int *l = link[j];
				float *p = vertex[j];
				for(int k=0; k<deg; k++){
					float *p1 = vertex[l[k]];
					float *p2 = vertex[l[(k+1)%deg]];
					
					float v1[3], v2[3];
					MeshData::VEC(v1, p, p1);
					MeshData::VEC(v2, p, p2);
					double c[3];
					MeshData::CROSS(c, v1, v2);
					double len = MeshData::LENGTH(c);
					if((float)len == 0)
						continue;

					float v3[3];
					MeshData::VEC(v3, p1, p2);

					double w1 = -MeshData::DOT(v1, v3)/len;
					double w2 = MeshData::DOT(v2, v3)/len;
					total_w += w1 + w2;
			
					L[j][0] += w2*v1[0] + w1*v2[0];
					L[j][1] += w2*v1[1] + w1*v2[1];
					L[j][2] += w2*v1[2] + w1*v2[2];
				}
				if(total_w != 0){
					L[j][0] /= total_w;
					L[j][1] /= total_w;
					L[j][2] /= total_w;
				}
				else
					L[j][0] = L[j][1] = L[j][2] = 0;
			}
		}
		
		int j;  for(i = 0; j<vertex_N; j++){
			float *p = mesh->vertex[j];
			//replace each vertex
			p[0] += (float)(dt*L[j][0]);
			p[1] += (float)(dt*L[j][1]);
			p[2] += (float)(dt*L[j][2]);
		}
	}
	delete[] L;
}

//void Smoother::moveNormalAreaD(float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal_f)[3] = mesh->normal_f;
//	double (*dp)[3] = new double[vertex_N][3];
//	int **link_f = mesh->vertex_link_f;
//	int **link_v = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//
//	int i;  for(i = 0; i<vertex_N; i++){
//		if(mesh->isBound[i]){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//
//		int deg = degree[i];
//		if(deg == 0){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//		int *l_f = link_f[i];
//		int *l_v = link_v[i];
//		
//		double df[3];
//		df[0] = df[1] = df[2] = 0;
//		float* p =vertex[i];
//		double total_w = 0;
//		for(int j = 0; j<deg; j++){
//			int f = l_f[j];
//			float* p1 = vertex[l_v[j]];
//			float* p2 = vertex[l_v[(j+1)%deg]];
//
//			float v1[3], v2[3], v3[3];
//			MeshData::VEC(v1, p, p1);
//			MeshData::VEC(v2, p, p2);
//			MeshData::VEC(v3, p1, p2);
//
//			double c[3];
//			MeshData::CROSS(c, v1, v2);
//			double A = MeshData::LENGTH(c);
//			if((float)A != 0){
//				double cot1 = -MeshData::DOT(v1, v3)/A;
//				double cot2 =  MeshData::DOT(v2, v3)/A;
//
//				df[0] += cot2*v1[0] + cot1*v2[0];
//				df[1] += cot2*v1[1] + cot1*v2[1];
//				df[2] += cot2*v1[2] + cot1*v2[2];
//			}
//
//			double m[3];
//			m[0] = normal_f[f][0];
//			m[1] = normal_f[f][1];
//			m[2] = normal_f[f][2];
//
//			double dot1 = MeshData::DOT(v1, m);
//			double dot2 = MeshData::DOT(v2, m);
//
//			double vp1[3], vp2[3], vp3[3];
//			vp1[0] = v1[0] - dot1*m[0];
//			vp1[1] = v1[1] - dot1*m[1];
//			vp1[2] = v1[2] - dot1*m[2];
//
//			vp2[0] = v2[0] - dot2*m[0];
//			vp2[1] = v2[1] - dot2*m[1];
//			vp2[2] = v2[2] - dot2*m[2];
//
//			vp3[0] = vp2[0] - vp1[0];
//			vp3[1] = vp2[1] - vp1[1];
//			vp3[2] = vp2[2] - vp1[2];
//
//			MeshData::CROSS(c, vp1, vp2);
//			A = MeshData::LENGTH(c);
//			if((float)A != 0){
//				double cot1 = -MeshData::DOT(vp1, vp3)/A;
//				double cot2 =  MeshData::DOT(vp2, vp3)/A;
//
//				df[0] -= cot2*vp1[0] + cot1*vp2[0];
//				df[1] -= cot2*vp1[1] + cot1*vp2[1];
//				df[2] -= cot2*vp1[2] + cot1*vp2[2];
//			}
//			
//			double w = 1;
//			total_w += w;
//		}
//		if(total_w != 0){
//			dp[i][0] = df[0]/total_w;
//			dp[i][1] = df[1]/total_w;
//			dp[i][2] = df[2]/total_w;
//		}
//		else
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//	}
//
//	for(i=0; i<vertex_N; i++){
//		vertex[i][0] += dt*(float)dp[i][0];
//		vertex[i][1] += dt*(float)dp[i][1];
//		vertex[i][2] += dt*(float)dp[i][2];
//	}
//	delete[] dp;
//}

void Smoother::smoothNormalV2(float T, float rate, int f_type)
{
	float (*vertex)[3] = mesh->vertex;
	float (*normal)[3] = mesh->normal_f;
	int (*face)[3] = mesh->face;
	int face_N = mesh->face_N;
	float (*tmp_normal)[3] = new float[face_N][3];
	int (*ad_face)[3] = mesh->face_link_E;
	int **ad_face2 = mesh->face_link_V;
	int *ad_face2_N = mesh->face_link_V_N;

	int i;  for(i = 0; i<face_N; i++){
		double total_w = 0;
		double sum_n[3];
		sum_n[0] = sum_n[1] = sum_n[2] = 0;
		double a = mesh->faceArea(i);
		/*
		double c[3];
		c[0] = (vertex[face[i][0]][0] + vertex[face[i][1]][0] + vertex[face[i][2]][0])/3.0;
		c[1] = (vertex[face[i][0]][1] + vertex[face[i][1]][1] + vertex[face[i][2]][1])/3.0;
		c[2] = (vertex[face[i][0]][2] + vertex[face[i][1]][2] + vertex[face[i][2]][2])/3.0;
		*/
		int j;  for(i = 0; j<3; j++){
			int pair = ad_face[i][j];
			if(pair < 0)
				continue;

			double al = mesh->faceArea(pair);
			if(al + a == 0)
				continue;

			double edge = MeshData::DIST(vertex[face[i][j]], vertex[face[i][(j+1)%3]]);
		
			double dot = MeshData::DOT(normal[i], normal[pair]);
			if(dot > 1)
				dot = 1;
			else if(dot < -1)
				dot = -1;

			double k = acos(dot)*edge/(a+al);;
			double w;
			if(f_type == 0)
				w = al/(1.0+(k/T)*(k/T));
			else
				w = al*exp(-(k/T)*(k/T));
			total_w += w;
			sum_n[0] += w*normal[pair][0];
			sum_n[1] += w*normal[pair][1];
			sum_n[2] += w*normal[pair][2];
		}
		/*
		for(j=0; j<ad_face2_N[i]; j++){
			int pair = ad_face2[i][j];
			if(pair < 0 || pair >= face_N)
				continue;
		
			double c1[3];
			c1[0] = (vertex[face[pair][0]][0] + vertex[face[pair][1]][0] + vertex[face[pair][2]][0])/3.0;
			c1[1] = (vertex[face[pair][0]][1] + vertex[face[pair][1]][1] + vertex[face[pair][2]][1])/3.0;
			c1[2] = (vertex[face[pair][0]][2] + vertex[face[pair][1]][2] + vertex[face[pair][2]][2])/3.0;
			double d = MeshData::DIST(c,c1);
			if((float)d == 0)
				continue;

			double area = MeshData::AREA(vertex[face[pair][0]], vertex[face[pair][1]], vertex[face[pair][2]]);
			
			double dot = MeshData::DOT(normal[i], normal[pair]);
			if(dot > 1)
				dot = 1;
			else if(dot < -1)
				dot = -1;

			double k = acos(dot)/d;

			double w;
			if(f_type == 0)
				w = rate*area/(1.0+(k/T)*(k/T)); 
			else
				w = rate*area*exp(-(k/T)*(k/T));
			total_w += w;
			sum_n[0] += w*normal[pair][0];
			sum_n[1] += w*normal[pair][1];
			sum_n[2] += w*normal[pair][2];
		}*/
		if((float)total_w != 0){
			tmp_normal[i][0] = (float)(sum_n[0]/total_w);
			tmp_normal[i][1] = (float)(sum_n[1]/total_w);
			tmp_normal[i][2] = (float)(sum_n[2]/total_w);

			double len = MeshData::LENGTH(tmp_normal[i]);
			if((float)len != 0){
				tmp_normal[i][0] = (float)(tmp_normal[i][0]/len);
				tmp_normal[i][1] = (float)(tmp_normal[i][1]/len);
				tmp_normal[i][2] = (float)(tmp_normal[i][2]/len);
			}
			else{
				tmp_normal[i][0] = normal[i][0];
				tmp_normal[i][1] = normal[i][1];
				tmp_normal[i][2] = normal[i][2];
			}
		}
		else{
			tmp_normal[i][0] = normal[i][0];
			tmp_normal[i][1] = normal[i][1];
			tmp_normal[i][2] = normal[i][2];
		}
	}
	delete[] normal;
	mesh->normal_f = tmp_normal;
}

//void Smoother::moveNormalAreaD(float dt, int power)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal_f)[3] = mesh->normal_f;
//	double (*dp)[3] = new double[vertex_N][3];
//	int **link_f = mesh->vertex_link_f;
//	int **link_v = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//
//	int i;  for(i = 0; i<vertex_N; i++){
//		if(mesh->isBound[i]){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//
//		int deg = degree[i];
//		if(deg == 0){
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//			continue;
//		}
//		int *l_f = link_f[i];
//		int *l_v = link_v[i];
//		
//		double df[3];
//		df[0] = df[1] = df[2] = 0;
//		float* p =vertex[i];
//		double total_w = 0;
//		for(int j = 0; j<deg; j++){
//			double grdA[3], grdB[3];
//			int f = l_f[j];
//			float* p1 = vertex[l_v[j]];
//			float* p2 = vertex[l_v[(j+1)%deg]];
//
//			float v1[3], v2[3], v3[3];
//			MeshData::VEC(v1, p, p1);
//			MeshData::VEC(v2, p, p2);
//			MeshData::VEC(v3, p1, p2);
//
//			double c[3];
//			MeshData::CROSS(c, v1, v2);
//			double A = MeshData::LENGTH(c);
//			if((float)A != 0){
//				double cot1 = -MeshData::DOT(v1, v3)/A;
//				double cot2 =  MeshData::DOT(v2, v3)/A;
//
//				grdA[0] = cot2*v1[0] + cot1*v2[0];
//				grdA[1] = cot2*v1[1] + cot1*v2[1];
//				grdA[2] = cot2*v1[2] + cot1*v2[2];
//			}
//			else{
//				grdA[0] = grdA[1] = grdA[2] = 0;
//			}
//
//			double m[3];
//			m[0] = normal_f[f][0];
//			m[1] = normal_f[f][1];
//			m[2] = normal_f[f][2];
//
//			double dot1 = MeshData::DOT(v1, m);
//			double dot2 = MeshData::DOT(v2, m);
//
//			double vp1[3], vp2[3], vp3[3];
//			vp1[0] = v1[0] - dot1*m[0];
//			vp1[1] = v1[1] - dot1*m[1];
//			vp1[2] = v1[2] - dot1*m[2];
//
//			vp2[0] = v2[0] - dot2*m[0];
//			vp2[1] = v2[1] - dot2*m[1];
//			vp2[2] = v2[2] - dot2*m[2];
//
//			vp3[0] = vp2[0] - vp1[0];
//			vp3[1] = vp2[1] - vp1[1];
//			vp3[2] = vp2[2] - vp1[2];
//
//			MeshData::CROSS(c, vp1, vp2);
//			double B = MeshData::LENGTH(c);
//			if((float)B != 0){
//				double cot1 = -MeshData::DOT(vp1, vp3)/B;
//				double cot2 =  MeshData::DOT(vp2, vp3)/B;
//
//				grdB[0] = cot2*vp1[0] + cot1*vp2[0];
//				grdB[1] = cot2*vp1[1] + cot1*vp2[1];
//				grdB[2] = cot2*vp1[2] + cot1*vp2[2];
//			}
//			else{
//				grdB[0] = grdB[1] = grdB[2] = 0;
//			}
//			
//			double w = (double)power*pow((A-B)/A, power-1.0);
//			total_w += w;
//			df[0] += w*(grdA[0] - grdB[0]);
//			df[1] += w*(grdA[1] - grdB[1]);
//			df[2] += w*(grdA[2] - grdB[2]);
//		}
//		if(total_w != 0){
//			dp[i][0] = df[0]/total_w;
//			dp[i][1] = df[1]/total_w;
//			dp[i][2] = df[2]/total_w;
//		}
//		else
//			dp[i][0] = dp[i][1] = dp[i][2] = 0;
//	}
//
//	for(i=0; i<vertex_N; i++){
//		vertex[i][0] += dt*(float)dp[i][0];
//		vertex[i][1] += dt*(float)dp[i][1];
//		vertex[i][2] += dt*(float)dp[i][2];
//	}
//	delete[] dp;
//}

//void Smoother::smoothRidgeEdges(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal)[3];
//	int *dist = new int[vertex_N];
//	int max_dist = 3;
//
//	Node** ridge = mesh->ridge_edge;
//	Node** ravine = mesh->ravine_edge;
//	Node** tag = new Node*[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++){
//		tag[i] = new Node;
//		for(Node* current = ridge[i]; current->next!=NULL; current=current->next)
//			tag[i]->append(current->v, -1);
//		for(current = ravine[i]; current->next!=NULL; current=current->next)
//			tag[i]->append(current->v, -1);
//		if(tag[i]->next!=NULL)
//			dist[i] = 0;
//		else
//			dist[i] = 10000;
//	}
//
//	for(i=0; i<max_dist; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(dist[j] == i){
//				int *nei, nei_N;
//				mesh->get1Ring(j, nei, nei_N);
//				if(nei_N == 0)
//					continue;
//				for(int k=0; k<nei_N; k++)
//					dist[nei[k]] = min(dist[nei[k]], i+1);
//			}
//		}
//	}
//
//	int tag_N = 0;
//	int near_N = 0;
//	for(i=0; i<vertex_N; i++){
//		if(dist[i] == 0)
//			tag_N++;
//		else if(dist[i] > 0)
//			near_N++;
//	}
//	int* tag_index = new int[tag_N];
//	int* tag_deg = new int[tag_N];
//	tag_N = 0;
//	int* near_index = new int[near_N];
//	near_N = 0;
//	for(i=0; i<vertex_N; i++){
//		if(dist[i] == 0){
//			tag_index[tag_N] = i;
//			int c = 0;
//			for(Node* current=tag[i]; current->next != NULL; current=current->next)
//				c++;
//			tag_deg[tag_N] = c;
//			tag_N++;
//		}
//		else if(dist[i] > 0){
//			near_index[near_N++] = i;
//		}
//	}
//
//
//	double (*L)[3] = new double[vertex_N][3];
//	double (*BL)[3] = new double[vertex_N][3];
//	for(i=0; i<iter; i++){
//		//smooth tag curve 
//		mesh->computeFaceNormal();
//		mesh->computeNormal();
//		normal = mesh->normal;
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(dist[j] != 0)
//				continue;
//			L[j][0] = L[j][1] = L[j][2] = 0;
//			int m = 0;
//			for(Node* current=tag[j]; current->next!=NULL; current=current->next){
//				m++;
//				L[j][0] += vertex[current->v][0];
//				L[j][1] += vertex[current->v][1];
//				L[j][2] += vertex[current->v][2];
//			}
//			if(m > 1){
//				L[j][0] = L[j][0]/(float)m - vertex[j][0];
//				L[j][1] = L[j][1]/(float)m - vertex[j][1];
//				L[j][2] = L[j][2]/(float)m - vertex[j][2];
//			
//				double inner = MeshData::DOT(L[j], normal[j]);
//				L[j][0] = L[j][0] - inner*normal[j][0];
//				L[j][1] = L[j][1] - inner*normal[j][1];
//				L[j][2] = L[j][2] - inner*normal[j][2];
//			}
//			else{
//				L[j][0] = L[j][1] = L[j][2] = 0;
//			}
//		}
//		for(j=0; j<vertex_N; j++){
//			if(dist[j] == 0){
//				vertex[j][0] += (float)(dt*L[j][0]);
//				vertex[j][1] += (float)(dt*L[j][1]);
//				vertex[j][2] += (float)(dt*L[j][2]);
//			}
//		}
//
//		for(j=0; j<vertex_N; j++){
//			L[j][0] = L[j][1] = L[j][2] = 0;
//			if(dist[j] > max_dist)
//				continue;
//
//			int *nei, nei_N;
//			mesh->get1Ring(j, nei, nei_N);
//			if(nei_N == 0)
//				continue;
//
//			for(int k=0; k<nei_N; k++){
//				L[j][0] += vertex[nei[k]][0];
//				L[j][1] += vertex[nei[k]][1];
//				L[j][2] += vertex[nei[k]][2];
//			}
//		
//			L[j][0] = L[j][0]/(float)nei_N - vertex[j][0];
//			L[j][1] = L[j][1]/(float)nei_N - vertex[j][1];
//			L[j][2] = L[j][2]/(float)nei_N - vertex[j][2];
//		}
//
//		for(j=0; j<vertex_N; j++){
//			BL[j][0] = BL[j][1] = BL[j][2] = 0;
//			if(dist[j] == 0 || dist[j] >= max_dist)
//				continue;
//	
//			int *nei, nei_N;
//			mesh->get1Ring(j, nei, nei_N);
//			if(nei_N == 0)
//				continue;
//
//			for(int k=0; k<nei_N; k++){
//				BL[j][0] += L[nei[k]][0];
//				BL[j][1] += L[nei[k]][1];
//				BL[j][2] += L[nei[k]][2];
//			}
//		
//			BL[j][0] = BL[j][0]/(float)nei_N - L[j][0];
//			BL[j][1] = BL[j][1]/(float)nei_N - L[j][1];
//			BL[j][2] = BL[j][2]/(float)nei_N - L[j][2];
//			/*
//			int m = 0;
//			for(Node* current=tag[j]; current->next!=NULL; current=current->next){
//				m++;
//				BL[j][0] += L[current->v][0];
//				BL[j][1] += L[current->v][1];
//				BL[j][2] += L[current->v][2];
//			}
//			if(m!=0){
//				BL[j][0] = BL[j][0]/(float)m - L[j][0];
//				BL[j][1] = BL[j][1]/(float)m - L[j][1];
//				BL[j][2] = BL[j][2]/(float)m - L[j][2];
//			}
//			*/
//		}
//
//		for(j=0; j<vertex_N; j++){
//			if(!mesh->isBound[j]){
//				vertex[j][0] -= (float)(dt*BL[j][0]);
//				vertex[j][1] -= (float)(dt*BL[j][1]);
//				vertex[j][2] -= (float)(dt*BL[j][2]);
//			}
//		}
//	}
//	delete[] L;
//	delete[] BL;
//	delete[] dist;
//	for(i=0; i<vertex_N; i++)
//		delete tag[i];
//	delete tag;
//}

//void Smoother::connectTag(int iter, float dt, float angle)
//{
//	float (*normal)[3] = mesh->normal_f;
//	float (*vertex)[3] = mesh->vertex;
//	Node** ridge_edge = mesh->ridge_edge;
//}

//void Smoother::Bilaplacian(int iter, float dt)
//{
//	int vertex_N = mesh->vertex_N;
//	double (*L)[3] = new double[vertex_N][3]; 
//	double (*B)[3] = new double[vertex_N][3]; 
//	BOOL *isBound = mesh->isBound;
//	int **link = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//	float (*vertex)[3] = mesh->vertex;
//
//	int i;  for(i = 0; i<iter; i++){
//		int j;  for(i = 0; j<vertex_N; j++){
//			if(isBound[j])
//				L[j][0] = L[j][1] = L[j][2] = 0;
//			else
//				mesh->laplacian(j, L[j]);
//		}
//		
//		for(j=0; j<vertex_N; j++){
//			B[j][0] = B[j][1] = B[j][2] = 0;
//			if(isBound[j])
//				continue;
//
//			int *l = link[j];
//			int deg = degree[j];
//			double v = 0;
//			for(int k=0; k<deg; k++){
//				int pair = l[k];
//				if(isBound[pair]){
//					v = 0;
//					break;
//				}
//				B[j][0] += L[j][0] - L[pair][0];
//				B[j][1] += L[j][1] - L[pair][1];
//				B[j][2] += L[j][2] - L[pair][2];
//				
//				v += 1.0/degree[pair];
//			}
//			if((float)v != 0){
//				v += deg;
//				B[j][0] /= v;
//				B[j][1] /= v;
//				B[j][2] /= v;
//			}
//			else
//				B[j][0] = B[j][1] = B[j][2] = 0;
//		}
//
//		for(j=0; j<vertex_N; j++){
//			vertex[j][0] += (float)(dt*B[j][0]);
//			vertex[j][1] += (float)(dt*B[j][1]);
//			vertex[j][2] += (float)(dt*B[j][2]);
//		}
//	}
//	
//	delete[] L;
//	delete[] B;
//}

void Smoother::smoothNormalMori(int iter)
{
	float (*normal)[3] = mesh->normal_f;
	int (*face)[3] = mesh->face;
	int face_N = mesh->face_N;
	float (*tmp_normal)[3] = new float[face_N][3];
	int (*ad_face)[3] = mesh->face_link_E;
	int count = 0;

	for(int k=0; k<iter; k++){
		count = 0;
		int i;  for(i = 0; i<face_N; i++){
			double sum_n[3];
			sum_n[0] = normal[i][0];
			sum_n[1] = normal[i][1];
			sum_n[2] = normal[i][2];

			float *n = normal[i];

			int *l = ad_face[i];
			float *m[3];
			m[0] = normal[l[0]];
			m[1] = normal[l[1]];
			m[2] = normal[l[2]];

			BOOL flag = true;
			int j;  for(i = 0; j<3; j++){
				double v[3];
				MeshData::CROSS(v, m[j], m[(j+1)%3]);
				double d1 = v[0]*(m[j][0] - n[0]) 
					+ v[1]*(m[j][1] - n[1]) + v[2]*(m[j][2] - n[2]);
				double d2 = v[0]*(m[j][0] - m[(j+2)%3][0]) 
					+ v[1]*(m[j][1] - m[(j+2)%3][1]) + v[2]*(m[j][2] - m[(j+2)%3][2]);
				if(d1*d2 < 0){
					flag = false;
					break;
				}
			}

			if(!flag){
				sum_n[0] = 0;
				sum_n[1] = 0;
				sum_n[2] = 0;

				for(j=0; j<3; j++){
					int pair = ad_face[i][j];
					if(pair < 0)
						continue;

					sum_n[0] += normal[pair][0];
					sum_n[1] += normal[pair][1];
					sum_n[2] += normal[pair][2];
				}
			}
			else{
				count++;
			}
			double len = MeshData::LENGTH(sum_n);
			if(len != 0){
				tmp_normal[i][0] = (float)(sum_n[0]/len);
				tmp_normal[i][1] = (float)(sum_n[1]/len);
				tmp_normal[i][2] = (float)(sum_n[2]/len);
			}
			else{
				tmp_normal[i][0] = normal[i][0];
				tmp_normal[i][1] = normal[i][1];
				tmp_normal[i][2] = normal[i][2];
			}
		}
		TRACE("%d\n", count);
		for(i=0; i<face_N; i++){
			normal[i][0] = tmp_normal[i][0];
			normal[i][1] = tmp_normal[i][1];
			normal[i][2] = tmp_normal[i][2];
		}
	}
	delete[] tmp_normal;
}

//void Smoother::badTriangles(BOOL *isBad)
//{
//	mesh->computeFaceNormal();
//	mesh->computeNormal();
//	float (*normal)[3] = mesh->normal_f;
//	int (*face)[3] = mesh->face;
//	int face_N = mesh->face_N;
//	int (*ad_face)[3] = mesh->face_link_E;
//	float (*normal_v)[3] = mesh->normal;
//
//	int count = 0;
//	int i;  for(i = 0; i<face_N; i++){
//		float *n = normal[i];
//
//		int *l = ad_face[i];
//		float *m[3];
//		//m[0] = normal[l[0]];
//		//m[1] = normal[l[1]];
//		//m[2] = normal[l[2]];
//
//		m[0] = normal_v[face[i][0]];
//		m[1] = normal_v[face[i][1]];
//		m[2] = normal_v[face[i][2]];
//
//		BOOL flag = true;
//		int j;  for(i = 0; j<3; j++){
//			double v[3];
//			MeshData::CROSS(v, m[j], m[(j+1)%3]);
//			double d1 = MeshData::DOT(v, n);
//				//v[0]*(m[j][0] - n[0]) 
//				//+ v[1]*(m[j][1] - n[1]) + v[2]*(m[j][2] - n[2]);
//			double d2 = MeshData::DOT(v, m[(j+2)%3]);
//				//v[0]*(m[j][0] - m[(j+2)%3][0]) 
//				//+ v[1]*(m[j][1] - m[(j+2)%3][1]) + v[2]*(m[j][2] - m[(j+2)%3][2]);
//			if(d1*d2 < 0){
//				count++;
//				flag = false;
//				break;
//			}
//		}
//		isBad[i] = !flag;
//	}
//	//TRACE("%d\n", count);
//}

void Smoother::smoothNormalMedian(BOOL isWeighted, int d_type)
{
	float (*normal)[3] = mesh->normal_f;
	int (*face)[3] = mesh->face;
	int face_N = mesh->face_N;
	float (*tmp_normal)[3] = new float[face_N][3];
	int (*ad_face)[3] = mesh->face_link_E;
	int **ad_face2 = mesh->face_link_V;
	int *ad_face2_N = mesh->face_link_V_N;
	mesh->computeCenter();
	float (*center)[3] = mesh->center;

	int i;  for(i = 0; i<face_N; i++){
		//count adjacent triangles
		int n = 0;
		int j;  for(i = 0; j<3; j++){
			int pair = ad_face[i][j];
			if(pair < 0)
				continue;
			n++;
		}
		if(isWeighted)
			n *= 2;
		for(j=0; j<ad_face2_N[i]; j++){
			int pair = ad_face2[i][j];
			if(pair < 0 || pair >= face_N)
				continue;
			n++;
		}

		//compute diffrences
		int *f_index = new int[n];
		double *diff = new double[n];
		float *nP = normal[i];
		n=0;
		for(j=0; j<3; j++){
			int pair = ad_face[i][j];
			if(pair < 0)
				continue;
			f_index[n] = pair;
			double dot = MeshData::DOT(nP, normal[pair]);
			if(dot > 1)
				dot = 1;
			if(dot < -1)
				dot = -1;
			diff[n] = acos(dot);
			if(d_type == 1){
				double dist = MeshData::DIST(center[i], center[pair]);
				if(dist != 0)
					diff[n] /= dist;
				else
					diff[n] = 0;
			}
			n++;
			if(isWeighted){
				diff[n] = diff[n-1];
				f_index[n] = f_index[n-1];
				n++;
			}
		}
		
		for(j=0; j<ad_face2_N[i]; j++){
			int pair = ad_face2[i][j];
			if(pair < 0 || pair >= face_N)
				continue;
			f_index[n] = pair;
			double dot = MeshData::DOT(nP, normal[pair]);
			if(dot > 1)
				dot = 1;
			if(dot < -1)
				dot = -1;
			diff[n] = acos(dot);
			if(d_type == 1){
				double dist = MeshData::DIST(center[i], center[pair]);
				if(dist != 0)
					diff[n] /= dist;
				else
					diff[n] = 0;
			}
			n++;
		}

		int m = n/2;
		if(m%2 != 0)
			m += 1;
		for(j=0; j<m; j++){
			int min = j;
			for(int k=j+1; k<n; k++)
				if(diff[min] > diff[k])
					min = k;
			int tmp = f_index[j];
			f_index[j] = f_index[min];
			f_index[min] = tmp;
			double tmp_d = diff[j];
			diff[j] = diff[min];
			diff[min] = tmp_d;
		}
		m = f_index[m-1];
		tmp_normal[i][0] = normal[m][0];
		tmp_normal[i][1] = normal[m][1];
		tmp_normal[i][2] = normal[m][2];

		delete[] f_index;
		delete[] diff;
	}
	delete[] normal;
	mesh->normal_f = tmp_normal;
}

//void Smoother::adaptiveSmoothingM(int iter, BOOL isWeighted, int int_type, int inner_iter, float int_step, int d_type, int ave_iter, int L_deg)
//{
//	int i,j;
//	
//	if(mesh->face_link_V == NULL)
//		mesh->generateFaceLinkV();
//
//	for(i=0; i<iter; i++){
//		mesh->computeFaceNormal();
//		
//		for(j=0; j<ave_iter; j++)
//			this->smoothNormalMedian(isWeighted, d_type);
//					
//		if(int_type == 0){
//			for(j=0; j<inner_iter; j++)
//				this->moveNormal(int_step);
//		}
//		else{
//			if(L_deg != 1){
//				for(j=0; j<inner_iter; j++)
//					this->moveNormalAreaD(int_step, L_deg);
//			}
//			else{
//				for(j=0; j<inner_iter; j++)
//					this->moveNormalAreaD(int_step);
//			}		
//		}
//	}
//}

//void Smoother::LaplacianFlowImplicit(float dt)
//{
//	PBCG *pbcg = new PBCG;
//	setupPBCG(pbcg);
//
//	unsigned long k;
//	int n = mesh->vertex_N;
//	int *degree = mesh->degree_v;
//	int **link = mesh->vertex_link_v;
//
//	double *sa = pbcg->sa;
//	unsigned long *ija = pbcg->ija;
//
//	for (int j=1;j<=n;j++) sa[j]= 1.0 + dt;
//	ija[1]=n+2;
//	k=n+1;
//	for(int i=1;i<=n;i++){
//		int deg = degree[i-1];
//		if(deg == 0 || mesh->isBound[i-1]){
//			sa[i] = 1.0;
//			ija[i+1]=k+1;
//			continue;
//		}
//		int *tmp = new int[deg];
//		int *l = link[i-1];
//		int j;  for(i = 0; j<deg; j++)
//			tmp[j] = l[j]+1;
//
//		//sorting indeces.
//		for(int s=0; s<deg; s++){
//			for(int t=s; t<deg; t++){
//				if(tmp[s] > tmp[t]){
//					int i_tmp = tmp[s];
//					tmp[s] = tmp[t];
//					tmp[t] = i_tmp;
//				}
//			}
//		}
//
//		for(j=0; j<deg; j++){
//			int index = tmp[j];
//			k++;
//			sa[k] = -dt/(double)deg;
//			ija[k] = index;
//		}
//		ija[i+1]=k+1;
//		delete[] tmp;
//	}
//
//	int iter;
//	double err;
//	double *old_vertex = new double[n+1];
//	double *new_vertex = new double[n+1];
//
//	float (*vertex)[3] = mesh->vertex;
//
//	//for x
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][0];
//		new_vertex[i+1] = vertex[i][0];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-3, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][0] = (float)new_vertex[i+1];
//
//	//for y
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][1];
//		new_vertex[i+1] = vertex[i][1];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-3, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][1] = (float)new_vertex[i+1];
//
//	//for z
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][2];
//		new_vertex[i+1] = vertex[i][2];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-3, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][2] = (float)new_vertex[i+1];
//
//	delete[] sa;
//	delete[] ija;
//	delete[] old_vertex;
//	delete[] new_vertex;
//	delete pbcg;
//}

//void Smoother::setupPBCG(PBCG *pbcg)
//{
//	int vertex_N = mesh->vertex_N;
//	int *degree = mesh->degree_v;
//	int size = vertex_N;
//	int i;  for(i = 0; i<vertex_N; i++)
//		size += degree[i];
//	
//	pbcg->sa = new double[size+2];
//	pbcg->ija = new unsigned long[size+2];
//}

//void Smoother::MeanCurvatureFlowImplicit(float dt)
//{
//	PBCG *pbcg = new PBCG;
//	setupPBCG(pbcg);
//
//	unsigned long k;
//	int n = mesh->vertex_N;
//	int *degree = mesh->degree_v;
//	int **link = mesh->vertex_link_v;
//	float (*vertex)[3] = mesh->vertex;
//	BOOL *isBound = mesh->isBound;
//
//	double *sa = pbcg->sa;
//	unsigned long *ija = pbcg->ija;
//
//	ija[1]=n+2;
//	k=n+1;
//	for(int i=1;i<=n;i++){
//		int ii = i-1;
//		int deg = degree[ii];
//		if(deg == 0 || isBound[ii]){
//			sa[i] = 1.0;
//			ija[i+1]=k+1;
//			continue;
//		}
//		int *l = link[ii];
//		double *w = new double[deg];
//		int j;  for(i = 0; j<deg; j++)
//			w[j] = 0;
//		double total_cot = 0;
//		double A = 0;
//		for(j=0; j<deg; j++){
//			int t = l[j];
//			int s = l[(j+1)%deg];
//
//			//vectors of edges of triangles
//			float PQ1[3], PQ2[3], QQ[3];
//			MeshData::VEC(PQ1, vertex[ii], vertex[s]);
//			MeshData::VEC(PQ2, vertex[ii], vertex[t]);
//			MeshData::VEC(QQ, vertex[s], vertex[t]);
//
//			//area of triangle
//			double Ai = MeshData::AREA(vertex[ii], vertex[s], vertex[t]);
//		
//			if(Ai > 0){ 
//				A += Ai;
//
//				double dot1 = -MeshData::DOT(PQ1, QQ);
//				double dot2 = MeshData::DOT(PQ2,QQ);
//
//				double cot1 = dot1/Ai;
//				double cot2 = dot2/Ai;
//
//				w[j] += cot1;
//				w[(j+1)%deg] += cot2;
//				total_cot += cot1 + cot2;
//			}
//		}
//		A *= 2.0;
//		if(A > 0){
//			sa[i] = 1.0 + dt*total_cot/A;
//			for(int k=0; k<deg; k++)
//				w[k] /= A;
//		}
//		else
//			sa[i] = 1.0;
//
//		int *tmp = new int[deg];
//		for(j=0; j<deg; j++)
//			tmp[j] = l[j]+1;
//
//		//sorting indeces.
//		for(int s=0; s<deg; s++){
//			for(int t=s; t<deg; t++){
//				if(tmp[s] > tmp[t]){
//					int i_tmp = tmp[s];
//					tmp[s] = tmp[t];
//					tmp[t] = i_tmp;
//					double w_tmp = w[s];
//					w[s] = w[t];
//					w[t] = w_tmp;
//				}
//			}
//		}
//
//		for(j=0; j<deg; j++){
//			if(w[j] == 0.0)
//				continue;
//			int index = tmp[j];
//			k++;
//			sa[k] = -dt*w[j];
//			ija[k] = index;
//		}
//		ija[i+1]=k+1;
//		delete[] tmp;
//		delete[] w;
//	}
//
//	int iter;
//	double err;
//	double *old_vertex = new double[n+1];
//	double *new_vertex = new double[n+1];
//
//	//for x
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][0];
//		new_vertex[i+1] = vertex[i][0];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-5, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][0] = (float)new_vertex[i+1];
//
//	//for y
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][1];
//		new_vertex[i+1] = vertex[i][1];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-5, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][1] = (float)new_vertex[i+1];
//
//	//for z
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][2];
//		new_vertex[i+1] = vertex[i][2];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-5, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][2] = (float)new_vertex[i+1];
//
//	delete[] old_vertex;
//	delete[] new_vertex;
//	delete pbcg;
//}

//void Smoother::IntegrateNormalImplicit(float dt)
//{
//	PBCG *pbcg = new PBCG;
//	setupPBCG(pbcg);
//
//	unsigned long k;
//	int n = mesh->vertex_N;
//	int *degree = mesh->degree_v;
//	int **link = mesh->vertex_link_v;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal)[3] = mesh->normal_f;
//	int **link_f = mesh->vertex_link_f;
//	BOOL *isBound = mesh->isBound;
//	float o[] = {0,0,0};
//	double *sa = pbcg->sa;
//	unsigned long *ija = pbcg->ija;
//
//	ija[1]=n+2;
//	k=n+1;
//	for(int i=1;i<=n;i++){
//		int ii = i-1;
//		int deg = degree[ii];
//		if(deg == 0 || isBound[ii]){
//			sa[i] = 1.0;
//			ija[i+1]=k+1;
//			continue;
//		}
//		int *l = link[ii];
//		int *l_f = link_f[ii];
//		double *w = new double[deg];
//		double *w_P = new double[deg];
//		int j;  for(i = 0; j<deg; j++){
//			w[j] = 0;
//			w_P[j] = 0;
//		}
//		double total_cot = 0;
//		double total_cot_P = 0;
//		for(j=0; j<deg; j++){
//			int t = l[j];
//			int s = l[(j+1)%deg];
//
//			//vectors of edges of triangles
//			float PQ1[3], PQ2[3], QQ[3];
//			MeshData::VEC(PQ1, vertex[ii], vertex[s]);
//			MeshData::VEC(PQ2, vertex[ii], vertex[t]);
//			MeshData::VEC(QQ, vertex[s], vertex[t]);
//
//			//area of triangle
//			double Ai = MeshData::AREA(vertex[ii], vertex[s], vertex[t]);
//		
//			if((float)Ai != 0){ 
//				double dot1 = -MeshData::DOT(PQ1, QQ);
//				double dot2 = MeshData::DOT(PQ2,QQ);
//
//				double cot1 = dot1/Ai;
//				double cot2 = dot2/Ai;
//
//				w[j] += cot1;
//				w[(j+1)%deg] += cot2;
//				total_cot += cot1 + cot2;
//			}
//
//			float* ni = normal[l_f[j]];
//			double dot = MeshData::DOT(ni, PQ1);
//			PQ1[0] = PQ1[0] - (float)dot*ni[0];
//			PQ1[1] = PQ1[1] - (float)dot*ni[1];
//			PQ1[2] = PQ1[2] - (float)dot*ni[2];
//
//			dot = MeshData::DOT(ni, PQ2);
//			PQ2[0] = PQ2[0] - (float)dot*ni[0];
//			PQ2[1] = PQ2[1] - (float)dot*ni[1];
//			PQ2[2] = PQ2[2] - (float)dot*ni[2];
//
//			MeshData::VEC(QQ, PQ1, PQ2);
//
//			Ai = MeshData::AREA(o, PQ1, PQ2);
//		
//			if((float)Ai != 0){ 
//				double dot1 = -MeshData::DOT(PQ1, QQ);
//				double dot2 = MeshData::DOT(PQ2,QQ);
//
//				double cot1 = dot1/Ai;
//				double cot2 = dot2/Ai;
//
//				w_P[j] += cot1;
//				w_P[(j+1)%deg] += cot2;
//				total_cot_P += cot1 + cot2;
//			}
//		}
//		
//		sa[i] = 1.0 + dt*(total_cot - total_cot_P);
//
//		int *tmp = new int[deg];
//		for(j=0; j<deg; j++)
//			tmp[j] = l[j]+1;
//
//		//sorting indeces.
//		for(int s=0; s<deg; s++){
//			for(int t=s; t<deg; t++){
//				if(tmp[s] > tmp[t]){
//					int i_tmp = tmp[s];
//					tmp[s] = tmp[t];
//					tmp[t] = i_tmp;
//					double w_tmp = w[s];
//					w[s] = w[t];
//					w[t] = w_tmp;
//					w_tmp = w_P[s];
//					w_P[s] = w_P[t];
//					w_P[t] = w_tmp;
//				}
//			}
//		}
//
//		for(j=0; j<deg; j++){
//			if(w[j] == 0 && w_P[j] == 0)
//				continue;
//			int index = tmp[j];
//			k++;
//			sa[k] = -dt*(w[j] - w_P[j]);
//			ija[k] = index;
//		}
//		ija[i+1]=k+1;
//		delete[] tmp;
//		delete[] w;
//	}
//
//	int iter;
//	double err;
//	double *old_vertex = new double[n+1];
//	double *new_vertex = new double[n+1];
//
//	//for x
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][0];
//		new_vertex[i+1] = vertex[i][0];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-3, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][0] = (float)new_vertex[i+1];
//
//	//for y
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][1];
//		new_vertex[i+1] = vertex[i][1];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-3, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][1] = (float)new_vertex[i+1];
//
//	//for z
//	for(i=0; i<n; i++){
//		old_vertex[i+1] = vertex[i][2];
//		new_vertex[i+1] = vertex[i][2];
//	}
//	pbcg->linbcg(n, old_vertex, new_vertex, 1, 1e-3, 100, &iter, &err);
//	for(i=0; i<n; i++)
//		vertex[i][2] = (float)new_vertex[i+1];
//
//	delete[] old_vertex;
//	delete[] new_vertex;
//	delete pbcg;
//}

//void Smoother::minimizeNormalErr()
//{
//	ThinPlate error;
//	//NormalError error;
//	error.mesh = mesh;
//	EnergyMinimizer mini;
//	mini.energy = &error;
//	
//	int n = mesh->vertex_N;
//	int iter1;
//	float fret;
//	float *p = new float[3*n];
//	int i;  for(i = 0; i<n; i++){
//		p[3*i] = mesh->vertex[i][0];
//		p[3*i+1] = mesh->vertex[i][1];
//		p[3*i+2] = mesh->vertex[i][2];
//	}
//
//	mini.frprmn(p, 3*n, 0.1f, &iter1, &fret);
//
//	for(i=0; i<n; i++){
//		mesh->vertex[i][0] = p[3*i];
//		mesh->vertex[i][1] = p[3*i+1];
//		mesh->vertex[i][2] = p[3*i+2];
//	}
//}

//void Smoother::smoothNormalGaussian(int size, float *sigma)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int *visit = new int[face_N];
//
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float d = (float)MeshData::DIST(center[i], center[ad]);
//				dist[i][j] = d;
//			}
//		}
//		visit[i] = -1;
//	}
//
//	for(i=0; i<face_N; i++){
//		double nor[] = {0,0,0};
//
//		Node* Q = new Node;
//		Q->append(i, 0);
//		visit[i] = i;
//		geo[i] = 0;
//		float *nf = normal[i];
//		while(Q->next != NULL){
//			int c = Q->v;
//			int r = Q->f;
//			Node* tmp = Q;
//			Q = Q->next;
//			Q->tail = tmp->tail;
//			tmp->next = NULL;
//			delete tmp;
//
//			if(geo[c] > 4*sqrt(sigma[i]))
//				continue;
//
//			//double diff = 1.5*(1.0-MeshData::DOT(normal[i], normal[c]));
//
//			double w;
//			if(geo[c] <= 2.0*sqrt(sigma[i]))
//				w = exp(-geo[c]*geo[c]/sigma[i]); // - diff);
//			else
//				w = 0.0625*pow(4.0 - geo[c]/sqrt(sigma[c]), 4.0)/(E*E);
//
//			nor[0] += area[c]*w*normal[c][0];
//			nor[1] += area[c]*w*normal[c][1];
//			nor[2] += area[c]*w*normal[c][2];
//
//			//if(r == size)
//				//continue;
//		
//			int* l = link[c];
//			int j;  for(i = 0; j<3; j++){
//				int k = l[j];
//				if(k < 0)
//					continue;
//				if(visit[k] == i)
//					continue;
//				Q->append(k, r+1);
//				visit[k] = i;
//
//				int* l_ad = link[k];
//				geo[k] = 100000;
//				double w1 = (1.0 - MeshData::DOT(nf,normal[k]));
//				//double b = 1.0+20.0*sqrt(1.0-dot);
//				float w = 5.0;
//				if(visit[l_ad[0]] == i){
//					double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[0]]));
//					double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[0]]))/3.0;
//					geo[k] = geo[l_ad[0]] + (1.0+w*(w1+w2+w3))*dist[k][0];
//				}
//				if(visit[l_ad[1]] == i){
//					double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[1]]));
//					double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[1]]))/3.0;
//					double d_tmp = geo[l_ad[1]] + (1.0+w*(w1+w2+w3))*dist[k][1];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//				if(visit[l_ad[2]] == i){
//					double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[2]]));
//					double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[2]]))/3.0;
//					double d_tmp = geo[l_ad[2]] + (1.0+w*(w1+w2+w3))*dist[k][2];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//
//
//				/*
//				float dot = MeshData::DOT(nf,normal[k]);
//				if(dot > 1.0)
//					dot = 1;
//				else if(dot < -1)
//					dot = -1;
//				double b = 1.0+20.0*sqrt(1.0-dot);
//				if(visit[l_ad[0]] == i)
//					geo[k] = geo[l_ad[0]] + b*dist[k][0];
//				if(visit[l_ad[1]] == i){
//					double d_tmp = geo[l_ad[1]] + b*dist[k][1];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//				if(visit[l_ad[2]] == i){
//					double d_tmp = geo[l_ad[2]] + b*dist[k][2];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}*/
//				//geo[k] = geo[c] + dist[c][j];
//			}
//		}
//		delete Q;
//		tmp_normal[i][0] = nor[0];
//		tmp_normal[i][1] = nor[1];
//		tmp_normal[i][2] = nor[2];
//	}
//
//	for(i=0; i<face_N; i++){
//		double len = MeshData::LENGTH(tmp_normal[i]);
//		if((float)len != 0){
//			normal[i][0] = (float)(tmp_normal[i][0]/len);
//			normal[i][1] = (float)(tmp_normal[i][1]/len);
//			normal[i][2] = (float)(tmp_normal[i][2]/len);
//		}
//	}
//
//	delete[] tmp_normal;
//	delete[] area;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//}

//void Smoother::decideSigma(float *sigma, float c, int size)
//{
//	int face_N = mesh->face_N;
//	int (*link)[3] = mesh->face_link_E;
//	float (*normal)[3] = mesh->normal_f;
//	int *visit = new int[face_N];
//	int i;  for(i = 0; i<face_N; i++)
//		visit[i] = -1;
//	
//	for(i=0; i<face_N; i++){
//		double M[] = {0,0,0};
//		int n = 0;
//
//		Node* Q = new Node;
//		Q->append(i, 0);
//		visit[i] = i;
//		while(Q->next != NULL){
//			int c = Q->v;
//			int r = Q->f;
//			Node* tmp = Q;
//			Q = Q->next;
//			Q->tail = tmp->tail;
//			tmp->next = NULL;
//			delete tmp;
//
//			M[0] += normal[c][0];
//			M[1] += normal[c][1];
//			M[2] += normal[c][2];
//			n++;
//
//			if(r == size)
//				continue;
//
//			int* l = link[c];
//			int j;  for(i = 0; j<3; j++){
//				int k = l[j];
//				if(k < 0)
//					continue;
//				if(visit[k] == i)
//					continue;
//				Q->append(k, r+1);
//				visit[k] = i;
//			}
//		}
//		delete Q;
//		M[0] /= n;
//		M[1] /= n;
//		M[2] /= n;
//		sigma[i] = c*c; //*(1.0 - MeshData::DOT(M, M));
//	}
//}

//void Smoother::checkGaussianSupport(int f, float sigma, float T)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int(*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//	if(mesh->triangle_id == NULL)
//		mesh->triangle_id = new int[face_N];
//	int *id = mesh->triangle_id;
//
//
//	int i;  for(i = 0; i<face_N; i++){
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0){
//				//dist[i][j] = 100000000;
//				continue;
//			}
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//			}
//		}
//		id[i] = -1;
//	}
//
//	Node* Q = new Node;
//	Q->append(f, 0);
//	id[f] = 0;
//	geo[f] = 0;
//	float *nf = normal[f];
//int count = 0;
//	while(Q->next != NULL){
//count++;
//		int c = Q->v;
//		int r = Q->f;
//		Node* tmp = Q;
//		Q = Q->next;
//		Q->tail = tmp->tail;
//		tmp->next = NULL;
//		delete tmp;
//
//		if(geo[c] > 4*sigma){
//			id[c] = -1;
//			continue;
//		}
//
//		//if(geo[c] > 2.0*sigma)
//			//id[c] = 2;
//		id[c] = (int)(255 - 255*geo[c]/(4*sigma));
//
//
//		int* l = link[c];
//		int j;  for(i = 0; j<3; j++){
//			int k = l[j];
//			if(k < 0)
//				continue;
//			if(id[k] >= 0)
//				continue;
//			Q->append(k, r+1);
//			id[k] = 1;
//
//			int* l_ad = link[k];
//			geo[k] = 100000;
//			double w1 = (1.0 - MeshData::DOT(nf,normal[k]));
//			//double b = 1.0+20.0*sqrt(1.0-dot);
//			if(l_ad[0] >= 0 && id[l_ad[0]] >= 0){
//				double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[0]]));
//				double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[0]]))/3.0;
//				geo[k] = geo[l_ad[0]] + (1.0+T*(w1+w2+w3))*dist[k][0];
//			}
//			if(l_ad[1] >= 0 && id[l_ad[1]] >= 0){
//				double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[1]]));
//				double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[1]]))/3.0;
//				double d_tmp = geo[l_ad[1]] + (1.0+T*(w1+w2+w3))*dist[k][1];
//				if(geo[k] > d_tmp)
//					geo[k] = d_tmp;
//			}
//			if(l_ad[2] >= 0 && id[l_ad[2]] >= 0){
//				double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[2]]));
//				double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[2]]))/3.0;
//				double d_tmp = geo[l_ad[2]] + (1.0+T*(w1+w2+w3))*dist[k][2];
//				if(geo[k] > d_tmp)
//					geo[k] = d_tmp;
//			}
//		}
//	}
//	delete Q;
//
//	delete[] dist;
//	delete[] geo;
//}

//float Smoother::computeSigma1(int f, int size)
//{
//	int face_N = mesh->face_N;
//	int (*link)[3] = mesh->face_link_E;
//	float (*normal)[3] = mesh->normal_f;
//	int *visit = new int[face_N];
//	int i;  for(i = 0; i<face_N; i++)
//		visit[i] = -1;
//	
//	double M[] = {0,0,0};
//	int n = 0;
//
//	Node* Q = new Node;
//	Q->append(f, 0);
//	visit[f] = f;
//	while(Q->next != NULL){
//		int c = Q->v;
//		int r = Q->f;
//		Node* tmp = Q;
//		Q = Q->next;
//		Q->tail = tmp->tail;
//		tmp->next = NULL;
//		delete tmp;
//
//		M[0] += normal[c][0];
//		M[1] += normal[c][1];
//		M[2] += normal[c][2];
//		n++;
//
//		if(r == size)
//			continue;
//
//		int* l = link[c];
//		int j;  for(i = 0; j<3; j++){
//			int k = l[j];
//			if(k < 0)
//				continue;
//			if(visit[k] == f)
//				continue;
//			Q->append(k, r+1);
//			visit[k] = f;
//		}
//	}
//	delete Q;
//	M[0] /= n;
//	M[1] /= n;
//	M[2] /= n;
//	return (float)(1.0 - MeshData::DOT(M, M));
//}

//float Smoother::computeSigma2(int f, int size)
//{
//	int face_N = mesh->face_N;
//	int (*link)[3] = mesh->face_link_E;
//	float (*normal)[3] = mesh->normal_f;
//	int *visit = new int[face_N];
//	int i;  for(i = 0; i<face_N; i++)
//		visit[i] = -1;
//	
//	double M = 0;
//	double sum = 0;
//	int n = 0;
//	float *nf = normal[f];
//
//	Node* Q = new Node;
//	Q->append(f, 0);
//	visit[f] = f;
//	while(Q->next != NULL){
//		int c = Q->v;
//		int r = Q->f;
//		Node* tmp = Q;
//		Q = Q->next;
//		Q->tail = tmp->tail;
//		tmp->next = NULL;
//		delete tmp;
//
//		double dot = MeshData::DOT(nf, normal[c]);
//		if(dot > 1.0)
//			dot = 1;
//		else if(dot < -1.0)
//			dot = -1.0;
//		double a = acos(dot);
//		M += a;
//		sum += a*a;
//		n++;
//
//		if(r == size)
//			continue;
//
//		int* l = link[c];
//		int j;  for(i = 0; j<3; j++){
//			int k = l[j];
//			if(k < 0)
//				continue;
//			if(visit[k] == f)
//				continue;
//			Q->append(k, r+1);
//			visit[k] = f;
//		}
//	}
//	delete Q;
//	M /= n;
//	sum /= n;
//
//	return (float)(sum - M*M);
//}

//void Smoother::computeOptimalGaussian(float s_min, float s_max, float step, float c)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int *visit = new int[face_N];
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//	//float c = 0.001;
//	//float c = 0.00005;
//
//	//float total_e = 0;
//	//int n = 0;
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//				//n++;
//				//total_e += dist[i][j];
//			}
//		}
//		visit[i] = -1;
//	}
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	for(i=0; i<face_N; i++){
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//		}
//
//		Node* Q = new Node;
//		Q->append(i, 0);
//		visit[i] = i;
//		geo[i] = 0;
//			
//		while(Q->next != NULL){
//			int c = Q->v;
//			int r = Q->f;
//			Node* tmp = Q;
//			Q = Q->next;
//			Q->tail = tmp->tail;
//			tmp->next = NULL;
//			delete tmp;
//
//			for(n=0; n<N; n++){
//				if(geo[c] > 4*si[n])
//					continue;
//				double w;
//				if(geo[c] <= 2.0*si[n])
//					w = area[c]*exp(-geo[c]*geo[c]/(si[n]*si[n])); 
//				else
//					w = area[c]*0.0625*pow(4.0 - geo[c]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[c][0];
//				ns[n][1] += w*normal[c][1];
//				ns[n][2] += w*normal[c][2];
//
//				total_w[n] += w;
//			}
//
//			if(geo[c] > 4*s_max)
//				continue;
//
//			int* l = link[c];
//			int j;  for(i = 0; j<3; j++){
//				int k = l[j];
//				if(k < 0)
//					continue;
//				if(visit[k] == i)
//					continue;
//				Q->append(k, r+1);
//				visit[k] = i;
//
//				int* l_ad = link[k];
//				geo[k] = 100000;
//				if(l_ad[0] >= 0 && visit[l_ad[0]] == i){
//					geo[k] = geo[l_ad[0]] + dist[k][0];
//				}
//				if(l_ad[1] >= 0 && visit[l_ad[1]] == i){
//					double d_tmp = geo[l_ad[1]] + dist[k][1];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//				if(l_ad[2] >= 0 && visit[l_ad[2]] == i){
//					double d_tmp = geo[l_ad[2]] + dist[k][2];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//			}
//		}
//		delete Q;
//
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			/*
//			double len = MeshData::LENGTH(ns[n]);
//			if((float)len == 0)
//				continue;
//			ns[n][0] /= len;
//			ns[n][1] /= len;
//			ns[n][2] /= len;
//			*/
//			//double v = c/(si[n]*si[n]) + 2.0*(1.0-MeshData::DOT(ni, ns[n]));
//			double v = c/(si[n]*si[n]) + 1.0-MeshData::LENGTH(ns[n])/total_w[n];
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns[opt]);
//		
//		tmp_normal[i][0] = ns[opt][0]/len;
//		tmp_normal[i][1] = ns[opt][1]/len;
//		tmp_normal[i][2] = ns[opt][2]/len;
//	}
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//
//	for(i=0; i<face_N; i++){
//		normal[i][0] = (float)tmp_normal[i][0];
//		normal[i][1] = (float)tmp_normal[i][1];
//		normal[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//	delete[] area;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//}

//bool Smoother::flipEdge(int i1, int i2, int &i3, int &i4)
//{
//	int* deg_v = mesh->degree_v;
//	int* deg_f = mesh->degree_f;
//	int **link_v = mesh->vertex_link_v;
//	int **link_f = mesh->vertex_link_f;
//	int (*link_e)[3] = mesh->face_link_E;
//	BOOL* isBound = mesh->isBound;
//	int (*face)[3] = mesh->face;
//
//	if(isBound[i1] || isBound[i2])
//		return false;
//
//	int deg1 = deg_v[i1];
//	int deg2 = deg_v[i2];
//
//	if(deg1 < 4 || deg2 < 4)
//		return false;
//
//	int* lv1 = link_v[i1];
//	int* lv2 = link_v[i2];
//	int* lf1 = link_f[i1];
//	int* lf2 = link_f[i2];
//
//	int f1, f2;
//	int i;  for(i = 0; i<deg1; i++){
//		if(lv1[i] == i2){
//			i3 = lv1[(i-1+deg1)%deg1];
//			i4 = lv1[(i+1)%deg1];
//			f1 = lf1[(i-1+deg1)%deg1];
//			f2 = lf1[i];
//			break;
//		}
//	}
//
//	if(isBound[i3] || isBound[i4])
//		return false;
//
//	int deg3 = deg_v[i3];
//	int deg4 = deg_v[i4];
//	int* lv3 = link_v[i3];
//	int* lv4 = link_v[i4];
//	int* lf3 = link_f[i3];
//	int* lf4 = link_f[i4];
//
//	for(i=0; i<deg3; i++){
//		if(lv3[i] == i4)
//			return false;
//	}
//
//	for(i=0; i<deg4; i++){
//		if(lv4[i] == i3)
//			return false;
//	}
//
//	//allocate new space 
//	int* n_lv1 = new int[deg1-1];
//	int* n_lv2 = new int[deg2-1];
//	int* n_lv3 = new int[deg3+1];
//	int* n_lv4 = new int[deg4+1];
//
//	int* n_lf1 = new int[deg1-1];
//	int* n_lf2 = new int[deg2-1];
//	int* n_lf3 = new int[deg3+1];
//	int* n_lf4 = new int[deg4+1];
//
//	//i1
//	int c = 0;
//	for(i=0; i<deg1; i++){
//		if(lv1[i] == i2)
//			continue;
//		n_lv1[c] = lv1[i];
//		n_lf1[c] = lf1[i];
//		c++;
//	}
//
//	//i2
//	c = 0;
//	for(i=0; i<deg2; i++){
//		if(lv2[i] == i1)
//			continue;
//		n_lv2[c] = lv2[i];
//		n_lf2[c] = lf2[i];
//		c++;
//	}
//
//	//i3
//	c = 0;
//	for(i=0; i<deg3; i++){
//		if(lv3[i] == i2){
//			n_lv3[c] = lv3[i];
//			n_lf3[c] = f2;
//			c++;
//
//			n_lv3[c] = i4;
//			n_lf3[c] = f1;
//			c++;
//		}
//		else{
//			n_lv3[c] = lv3[i];
//			n_lf3[c] = lf3[i];
//			c++;
//		}
//	}
//
//	//i4
//	c = 0;
//	for(i=0; i<deg4; i++){
//		if(lv4[i] == i1){
//			n_lv4[c] = lv4[i];
//			n_lf4[c] = f1;
//			c++;
//
//			n_lv4[c] = i3;
//			n_lf4[c] = f2;
//			c++;
//		}
//		else{
//			n_lv4[c] = lv4[i];
//			n_lf4[c] = lf4[i];
//			c++;
//		}
//	}
//
//	delete[] lv1;
//	delete[] lv2;
//	delete[] lv3;
//	delete[] lv4;
//
//	delete[] lf1;
//	delete[] lf2;
//	delete[] lf3;
//	delete[] lf4;
//
//	//replace link data
//	link_v[i1] = n_lv1;
//	link_v[i2] = n_lv2;
//	link_v[i3] = n_lv3;
//	link_v[i4] = n_lv4;
//
//	link_f[i1] = n_lf1;
//	link_f[i2] = n_lf2;
//	link_f[i3] = n_lf3;
//	link_f[i4] = n_lf4;
//
//	//update degree
//	deg_v[i1]--;
//	deg_v[i2]--;
//	deg_v[i3]++;
//	deg_v[i4]++;
//
//	deg_f[i1]--;
//	deg_f[i2]--;
//	deg_f[i3]++;
//	deg_f[i4]++;
//
//	//update face 
//	for(i=0; i<3; i++){
//		if(face[f1][i] == i2)
//			face[f1][i] = i4;
//	}
//
//	for(i=0; i<3; i++){
//		if(face[f2][i] == i1)
//			face[f2][i] = i3;
//	}
//
//	//update face normal
//	float (*normal)[3] = mesh->normal_f;
//	float (*vertex)[3] = mesh->vertex;
//	float v1[3], v2[3];
//	int *f = face[f1];
//	MeshData::VEC(v1, vertex[f[0]], vertex[f[1]]);
//	MeshData::VEC(v2, vertex[f[0]], vertex[f[2]]);
//	double n[3];
//	MeshData::CROSS(n, v1, v2);
//	double len = MeshData::LENGTH(n);
//	if((float)len != 0){
//		normal[f1][0] = (float)(n[0]/len);
//		normal[f1][1] = (float)(n[1]/len);
//		normal[f1][2] = (float)(n[2]/len);
//	}
//	else
//		normal[f1][0] = normal[f1][1] = normal[f1][2] = 0;
//
//	f = face[f2];
//	MeshData::VEC(v1, vertex[f[0]], vertex[f[1]]);
//	MeshData::VEC(v2, vertex[f[0]], vertex[f[2]]);
//	MeshData::CROSS(n, v1, v2);
//	len = MeshData::LENGTH(n);
//	if((float)len != 0){
//		normal[f2][0] = (float)(n[0]/len);
//		normal[f2][1] = (float)(n[1]/len);
//		normal[f2][2] = (float)(n[2]/len);
//	}
//	else
//		normal[f2][0] = normal[f2][1] = normal[f2][2] = 0;
//
//	//update face link
//	int *e1 = link_e[f1];
//	int *e2 = link_e[f2];
//	int index1, index2, f12, f21;
//	for(i=0; i<3; i++){
//		if(e1[i] == f2){
//			index1 = (i+2)%3;
//			f12 = e1[index1];
//			e1[index1] = f2;
//			index1 = i;
//			break;
//		}
//	}
//	for(i=0; i<3; i++){
//		if(e2[i] == f1){
//			index2 = (i+2)%3;
//			f21 = e2[index2];
//			e2[index2] = f1;
//			index2 = i;
//			break;
//		}
//	}
//	e1[index1] = f21;
//	e2[index2] = f12;
//
//	e1 = link_e[f12];
//	for(i=0; i<3; i++){
//		if(e1[i] == f1){
//			e1[i] = f2;
//			break;
//		}
//	}
//	e2 = link_e[f21];
//	for(i=0; i<3; i++){
//		if(e2[i] == f2){
//			e2[i] = f1;
//			break;
//		}
//	}
//
//	return true;
//}

//void Smoother::antialiasingEdgeFlip()
//{
//	int face_N = mesh->face_N;
//	int (*face)[3] = mesh->face;
//	int (*f_link)[3] = mesh->face_link_E;
//	float (*normal)[3] = mesh->normal_f;
//
//	int vertex_N = mesh->vertex_N;
//	int **e_table = new int*[vertex_N];
//	int *degree = mesh->degree_v;
//	int *edge_count = new int[vertex_N];
//	int i;  for(i = 0; i<vertex_N; i++){
//		e_table[i] = new int[degree[i]];
//		edge_count[i] = 0;
//	}
//
//	int edge_N = 0;
//	for(i=0; i<face_N; i++){
//		int *f = face[i];
//		int j;  for(i = 0; j<3; j++){
//			if(f_link[i][j] < 0)
//				continue;
//			edge_N++;
//		}
//	}
//	edge_N /= 2;
//
//	int (*edge)[2] = new int[edge_N][2];
//	float *q = new float[edge_N];
//	edge_N = 0;
//	for(i=0; i<face_N; i++){
//		int *f = face[i];
//		int j;  for(i = 0; j<3; j++){
//			if(f_link[i][j] < 0)
//				continue;
//			int i1 = f[j];
//			int i2 = f[(j+1)%3];
//			if(i1< i2){
//				edge[edge_N][0] = i1;
//				edge[edge_N][1] = i2;
//				q[edge_N] = -computeAliasingMeasure(i1, i2);
//				e_table[i1][edge_count[i1]++] = edge_N;
//				e_table[i2][edge_count[i2]++] = edge_N;
//				edge_N++;
//			}
//		}
//	}
//	delete[] edge_count;
//
//	//Construct priority Q
//	int* heap = new int[edge_N+1];
//	int* index = new int[edge_N];
//	int last_heap = 0;
//	for(i=0; i<edge_N; i++){
//		if(q[i] > -0.0000001){
//			index[i] = -1;
//			continue;
//		}
//		//insert;
//		heap[++last_heap] = i;
//		index[i] = last_heap;
//		upheap(q, last_heap, last_heap, heap, index);
//	}
//
//	while(last_heap != 0){
//		int min = heap[1];
//		index[min] = -1;
//		//remove
//		if(last_heap != 1){
//			heap[1] = heap[last_heap--];
//			index[heap[1]] = 1;
//			downheap(q, last_heap, 1, heap, index);
//		}
//		else
//			last_heap = 0;
//		int i1 = edge[min][0];
//		int i2 = edge[min][1];
//		int i3, i4;
//
//		if(!flipEdge(i1, i2, i3, i4))
//			continue;
//
//		edge[min][0] = i3;
//		edge[min][1] = i4;
//
//		int deg1 = degree[i1];
//		int *oe1 = e_table[i1];
//		int *ne1 = new int[deg1];
//		int c = 0;
//		int i;  for(i = 0; i<deg1+1; i++){
//			if(oe1[i] == min)
//				continue;
//			ne1[c++] = oe1[i];
//		}
//		delete[] oe1;
//		e_table[i1] = ne1;
//		for(i=0; i<deg1; i++){
//			int e = ne1[i];
//			int i1 = edge[e][0];
//			int i2 = edge[e][1];
//			float tmp = -computeAliasingMeasure(i1, i2);
//			if(q[e] == tmp)
//				continue;
//			q[e] = tmp;
//			if(index[e] < 0){
//				if(q[e] > -0.0000001)
//					continue;
//				else{
//					//continue;
//					//insert into Q
//					heap[++last_heap] = e;
//					index[e] = last_heap;
//					upheap(q, last_heap, last_heap, heap, index);
//				}
//			}
//			else{
//				if(q[e] > -0.0000001){
//					//delete from Q
//					float tmp = q[e];
//					q[e] = -1000000;
//					upheap(q, last_heap, index[e], heap, index);
//					index[e] = -1;
//					heap[1] = heap[last_heap--];
//					index[heap[1]] = 1;
//					downheap(q, last_heap, 1, heap, index);
//					q[e] = tmp;
//				}
//				else{
//					//update Q
//					if(index[e] != 1 && q[e] < q[heap[index[e]/2]])
//						upheap(q, last_heap, index[e], heap, index);
//					else
//						downheap(q, last_heap, index[e], heap, index);
//				}
//			}
//		}
//		
//		int deg2 = degree[i2];
//		int *oe2 = e_table[i2];
//		int *ne2 = new int[deg2];
//		c = 0;
//		for(i=0; i<deg2+1; i++){
//			if(oe2[i] == min)
//				continue;
//			ne2[c++] = oe2[i];
//		}
//		delete[] oe2;
//		e_table[i2] = ne2;
//		for(i=0; i<deg2; i++){
//			int e = ne2[i];
//			int i1 = edge[e][0];
//			int i2 = edge[e][1];
//			float tmp = -computeAliasingMeasure(i1, i2);
//			if(q[e] == tmp)
//				continue;
//			q[e] = tmp;
//			if(index[e] < 0){
//				if(q[e] > -0.0000001)
//					continue;
//				else{
//					//continue;
//					//insert into Q
//					heap[++last_heap] = e;
//					index[e] = last_heap;
//					upheap(q, last_heap, last_heap, heap, index);
//				}
//			}
//			else{
//				if(q[e] > -0.0000001){
//					//delete from Q
//					float tmp = q[e];
//					q[e] = -1000000;
//					upheap(q, last_heap, index[e], heap, index);
//					index[e] = -1;
//					heap[1] = heap[last_heap--];
//					index[heap[1]] = 1;
//					downheap(q, last_heap, 1, heap, index);
//					q[e] = tmp;
//				}
//				else{
//					//update Q
//					if(index[e] != 1 && q[e] < q[heap[index[e]/2]])
//						upheap(q, last_heap, index[e], heap, index);
//					else
//						downheap(q, last_heap, index[e], heap, index);
//				}
//			}
//		}
//
//		int deg3 = degree[i3];
//		int *oe3 = e_table[i3];
//		int *ne3 = new int[deg3];
//		for(i=0; i<deg3-1; i++)
//			ne3[i] = oe3[i];
//		ne3[i] = min;
//		delete[] oe3;
//		e_table[i3] = ne3;
//		for(i=0; i<deg3-1; i++){
//			int e = ne3[i];
//			int i1 = edge[e][0];
//			int i2 = edge[e][1];
//			float tmp = -computeAliasingMeasure(i1, i2);
//			if(q[e] == tmp)
//				continue;
//			q[e] = tmp;
//			if(index[e] < 0){
//				if(q[e] > -0.0000001)
//					continue;
//				else{
//					//continue;
//					//insert into Q
//					heap[++last_heap] = e;
//					index[e] = last_heap;
//					upheap(q, last_heap, last_heap, heap, index);
//				}
//			}
//			else{
//				if(q[e] > -0.0000001){
//					//delete from Q
//					float tmp = q[e];
//					q[e] = -1000000;
//					upheap(q, last_heap, index[e], heap, index);
//					index[e] = -1;
//					heap[1] = heap[last_heap--];
//					index[heap[1]] = 1;
//					downheap(q, last_heap, 1, heap, index);
//					q[e] = tmp;
//				}
//				else{
//					//update Q
//					if(index[e] != 1 && q[e] < q[heap[index[e]/2]])
//						upheap(q, last_heap, index[e], heap, index);
//					else
//						downheap(q, last_heap, index[e], heap, index);
//				}
//			}
//		}
//
//		int deg4 = degree[i4];
//		int *oe4 = e_table[i4];
//		int *ne4 = new int[deg4];
//		for(i=0; i<deg4-1; i++)
//			ne4[i] = oe4[i];
//		ne4[i] = min;
//		delete[] oe4;
//		e_table[i4] = ne4;
//		for(i=0; i<deg4-1; i++){
//			int e = ne4[i];
//			int i1 = edge[e][0];
//			int i2 = edge[e][1];
//			float tmp = -computeAliasingMeasure(i1, i2);
//			if(q[e] == tmp)
//				continue;
//			q[e] = tmp;
//			if(index[e] < 0){
//				if(q[e] > -0.0000001)
//					continue;
//				else{
//					//continue;
//					//insert into Q
//					heap[++last_heap] = e;
//					index[e] = last_heap;
//					upheap(q, last_heap, last_heap, heap, index);
//				}
//			}
//			else{
//				if(q[e] > -0.0000001){
//					//delete from Q
//					float tmp = q[e];
//					q[e] = -1000000;
//					upheap(q, last_heap, index[e], heap, index);
//					index[e] = -1;
//					heap[1] = heap[last_heap--];
//					index[heap[1]] = 1;
//					downheap(q, last_heap, 1, heap, index);
//					q[e] = tmp;
//				}
//				else{
//					//update Q
//					if(index[e] != 1 && q[e] < q[heap[index[e]/2]])
//						upheap(q, last_heap, index[e], heap, index);
//					else
//						downheap(q, last_heap, index[e], heap, index);
//				}
//			}
//		}
//	}
//
//	delete[] q;
//	delete[] edge;
//	delete[] heap;
//	delete[] index;
//	for(i=0; i<vertex_N; i++){
//		if(degree[i] != 0)
//			delete[] e_table[i];
//	}
//	delete[] e_table;
//}

inline void Smoother::upheap(float *a, int N, int k, int *p, int *q)
{
	int v;
	v = p[k];
	while(k > 1 && a[p[k/2]] >= a[v]){
		p[k] = p[k/2]; q[p[k/2]] = k; k = k/2;
	}
	p[k] = v; q[v] = k;
}

//inline void Smoother::downheap(float *a, int N, int k, int *p, int *q)
//{
//	int j, v;
//	v = p[k];
//	while(k <= N/2){
//		j = k+k;
//		if(j < N && a[p[j]] > a[p[j+1]]) j++;
//		if(a[v] <= a[p[j]]) break;
//		p[k] = p[j]; q[p[j]] = k; k = j;
//	}
//	p[k] = v; q[v] = k;
//}

//float Smoother::computeAliasingMeasure(int i1, int i2)
//{
//	int f1, index1;
//	int deg = mesh->degree_f[i1];
//	int *link_v = mesh->vertex_link_v[i1];
//	int *link_f = mesh->vertex_link_f[i1];
//	int i;  for(i = 0; i<deg; i++){
//		if(link_v[i] == i2){
//			f1 = link_f[i];
//			break;
//		}
//	}
//	int *v1 = mesh->face[f1];
//	if(v1[0] == i1){
//		if(v1[1] == i2)
//			index1 = 0;
//		else
//			index1 = 2;
//	}
//	else if(v1[1] == i1){
//		if(v1[0] == i2)
//			index1 = 0;
//		else
//			index1 = 1;
//	}
//	else{
//		if(v1[0] == i2)
//			index1 = 2;
//		else
//			index1 = 1;
//	}
//
//	int *l1 = mesh->face_link_E[f1];
//	int f2 = l1[index1];
//	if(f2 < 0)
//		return 0;
//	int *v2 = mesh->face[f2];
//	int *l2 = mesh->face_link_E[f2];
//	
//	float *p1 = mesh->vertex[v1[index1]];
//	float *p2 = mesh->vertex[v1[(index1+1)%3]];
//	float *p3 = mesh->vertex[v1[(index1+2)%3]];
//	float *p4;
//	int index2;
//	if(l2[0] == f1){
//		p4 = mesh->vertex[v2[2]];
//		index2 = 0;
//	}
//	else if(l2[1] == f1){
//		p4 = mesh->vertex[v2[0]];
//		index2 = 1;
//	}
//	else {
//		p4 = mesh->vertex[v2[1]];
//		index2 = 2;
//	}
//
//	float (*normal)[3] = mesh->normal_f;
//
//	float *n1 = normal[f1];
//	if(l1[(index1+1)%3] < 0)
//		return 0;
//	float *n1f = normal[l1[(index1+1)%3]];
//	if(l1[(index1+2)%3] < 0)
//		return 0;
//	float *n1b = normal[l1[(index1+2)%3]];
//
//	float *n2 = normal[f2];
//	if(l2[(index2+2)%3] < 0)
//		return 0;
//	float *n2f = normal[l2[(index2+2)%3]];
//	if(l2[(index2+1)%3] < 0)
//		return 0;
//	float *n2b = normal[l2[(index2+1)%3]];
//
//	float b1, b2, b1m, b2m;
//	b1 = b2 = b1m = b2m = (float)MeshData::DOT(n1, n2);
//	float dot = (float)MeshData::DOT(n1, n1f);
//	b1 += dot;
//	b1m = min(b1m, dot);
//	dot = (float)MeshData::DOT(n1, n1b);
//	b1 += dot;
//	b1m = min(b1m, dot);
//	//b1 -= b1m;
//
//	dot = (float)MeshData::DOT(n2, n2f);
//	b2 += dot;
//	b2m = min(b2m, dot);
//	dot = (float)MeshData::DOT(n2, n2b);
//	b2 += dot;
//	b2m = min(b2m, dot);
//	//b2 -= b2m;
//
//	float m1[3], m2[3];
//	float e1[3], e2[3];
//	MeshData::VEC(e1, p2, p3);
//	MeshData::VEC(e2, p2, p4);
//	double n[3];
//	MeshData::CROSS(n, e1, e2);
//	double len = MeshData::LENGTH(n);
//	if((float)len != 0){
//		m1[0] = (float)(n[0]/len);
//		m1[1] = (float)(n[1]/len);
//		m1[2] = (float)(n[2]/len);
//	}
//	else
//		return 0;
//	
//	if(MeshData::DOT(n1,m1) < 0 || MeshData::DOT(n2,m1) < 0)
//		return 0;
//
//	MeshData::VEC(e1, p1, p4);
//	MeshData::VEC(e2, p1, p3);
//	n[3];
//	MeshData::CROSS(n, e1, e2);
//	len = MeshData::LENGTH(n);
//	if((float)len != 0){
//		m2[0] = (float)(n[0]/len);
//		m2[1] = (float)(n[1]/len);
//		m2[2] = (float)(n[2]/len);
//	}
//	else
//		return 0;
//
//	if(MeshData::DOT(n1,m2) < 0 || MeshData::DOT(n2,m2) < 0)
//		return 0;
//
//	float a1, a2, a1m, a2m;
//	a1 = a2 = a1m = a2m = (float)MeshData::DOT(m1, m2);
//	dot = (float)MeshData::DOT(m1, n1f);
//	a1 += dot;
//	a1m = min(a1m, dot);
//	dot = (float)MeshData::DOT(m1, n2f);
//	a1 += dot;
//	a1m = min(a1m, dot);
//	//a1 -= a1m;
//
//	dot = (float)MeshData::DOT(m2, n1b);
//	a2 += dot;
//	a2m = min(a2m, dot);
//	dot = (float)MeshData::DOT(m2, n2b);
//	a2 += dot;
//	a2m = min(a2m, dot);
//	//a2 -= a2m;
//
//	return (a1 + a2) - (b1 + b2);
//}

//void Smoother::adaptiveGaussianSmoothing(float c, BOOL is_anistoropic, float w)
//{
//	float e = mesh->averageOfEdgeLength();
//
//	/***********************************
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	float (*o_normal)[3] = new float[face_N][3];;
//	int i;  for(i = 0; i<face_N; i++){
//		o_normal[i][0] = normal[i][0];
//		o_normal[i][1] = normal[i][1];
//		o_normal[i][2] = normal[i][2];
//	}
//	***********************************/
//	
//	if(!is_anistoropic)
//		//computeOptimalGaussian(0.25f*e, 2.5f*e, 0.25f*e, c*e*e);
//		computeOptimalGaussianDikstra(0.2f*e, 2.0f*e, 0.2f*e, c*e*e);
//	else
//		//computeAnistoropicGaussian(0.5f*e, 5.0f*e, 0.5f*e, 10.0f*c*e*e, w/e);
//		//computeAnistoropicGaussianDikstra(0.3f*e, 3.0f*e, 0.3f*e, 10.0f*c*e*e, w);
//		computeAnistoropicGaussianDikstra(0.4f*e, 4.0f*e, 0.4f*e, 10.0f*c*e*e, w);
//		//computeAnistoropicGaussianDikstra2(0.4f*e, 4.0f*e, 0.4f*e, 10*c*e*e, e*e*w);
//		//computeAnistoropicGaussianDikstraAngle(0.3f*e, 3.0f*e, 0.3f*e, 10.0f*c*e*e, w);
//
//	/***********************************
//	for(i=0; i<face_N; i++){
//		normal[i][0] = o_normal[i][0] + 2.5f*(o_normal[i][0] - normal[i][0]);
//		normal[i][1] = o_normal[i][1] + 2.5f*(o_normal[i][1] - normal[i][1]);
//		normal[i][2] = o_normal[i][2] + 2.5f*(o_normal[i][2] - normal[i][2]);
//		double len = MeshData::LENGTH(normal[i]);
//		if((float)len == 0)
//			continue;
//		normal[i][0] = (float)(normal[i][0]/len);
//		normal[i][1] = (float)(normal[i][1]/len);
//		normal[i][2] = (float)(normal[i][2]/len);
//	}
//	***********************************/
//	
//	minimizeNormalErr();
//}

//void Smoother::computeAnistoropicGaussian(float s_min, float s_max, float step, float c, float T)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int *visit = new int[face_N];
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//	//float c = 0.001;
//	//float c = 0.00005;
//
//	//float total_e = 0;
//	//int n = 0;
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//
//				//float d = (float)MeshData::DIST(center[i], center[ad]);
//				//dist[i][j] = d;
//				//n++;
//				//total_e += dist[i][j];
//			}
//		}
//		visit[i] = -1;
//	}
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	for(i=0; i<face_N; i++){
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//		}
//
//		Node* Q = new Node;
//		Q->append(i, 0);
//		visit[i] = i;
//		geo[i] = 0;
//		float *nf = normal[i];
//		
//		while(Q->next != NULL){
//			int c = Q->v;
//			int r = Q->f;
//			Node* tmp = Q;
//			Q = Q->next;
//			Q->tail = tmp->tail;
//			tmp->next = NULL;
//			delete tmp;
//
//			for(n=0; n<N; n++){
//				if(geo[c] > 4*si[n])
//					continue;
//				double w;
//				if(geo[c] <= 2.0*si[n])
//					w = area[c]*exp(-geo[c]*geo[c]/(si[n]*si[n])); 
//				else
//					w = area[c]*0.0625*pow(4.0 - geo[c]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[c][0];
//				ns[n][1] += w*normal[c][1];
//				ns[n][2] += w*normal[c][2];
//
//				total_w[n] += w;
//			}
//
//			if(geo[c] > 4*s_max)
//				continue;
//
//			int* l = link[c];
//			int j;  for(i = 0; j<3; j++){
//				int k = l[j];
//				if(k < 0)
//					continue;
//				if(visit[k] == i)
//					continue;
//				Q->append(k, r+1);
//				visit[k] = i;
//
//
//				int* l_ad = link[k];
//				geo[k] = 100000;
//				double w1 = (1.0 - MeshData::DOT(nf,normal[k]));
//				//double b = 1.0+20.0*sqrt(1.0-dot);
//				if(l_ad[0] >= 0 && visit[l_ad[0]] == i){
//					double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[0]]));
//					double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[0]]))/3.0;
//					geo[k] = geo[l_ad[0]] + (1.0+T*(w1+w2+w3))*dist[k][0];
//				}
//				if(l_ad[1] >= 0 && visit[l_ad[1]] == i){
//					double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[1]]));
//					double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[1]]))/3.0;
//					double d_tmp = geo[l_ad[1]] + (1.0+T*(w1+w2+w3))*dist[k][1];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//				if(l_ad[2] >= 0 && visit[l_ad[2]] == i){
//					double w2 = (1.0 - MeshData::DOT(nf,normal[l_ad[2]]));
//					double w3 = (1.0 - MeshData::DOT(normal[k],normal[l_ad[2]]))/3.0;
//					double d_tmp = geo[l_ad[2]] + (1.0+T*(w1+w2+w3))*dist[k][2];
//					if(geo[k] > d_tmp)
//						geo[k] = d_tmp;
//				}
//			}
//		}
//		delete Q;
//
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			/*
//			double len = MeshData::LENGTH(ns[n]);
//			if((float)len == 0)
//				continue;
//			ns[n][0] /= len;
//			ns[n][1] /= len;
//			ns[n][2] /= len;
//			*/
//			
//			//double v = c/(si[n]*si[n]) + 2.0*(1.0-MeshData::DOT(ni, ns[n]));
//			double v = c/(si[n]*si[n]) + 1.0-MeshData::LENGTH(ns[n])/total_w[n];
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns[opt]);
//		
//		tmp_normal[i][0] = ns[opt][0]/len;
//		tmp_normal[i][1] = ns[opt][1]/len;
//		tmp_normal[i][2] = ns[opt][2]/len;
//	}
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//
//	for(i=0; i<face_N; i++){
//		normal[i][0] = (float)tmp_normal[i][0];
//		normal[i][1] = (float)tmp_normal[i][1];
//		normal[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//	delete[] area;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//}

//void Smoother::addEdgeLengthNoise(float rate)
//{
//	int vertex_N = mesh->vertex_N;
//	float (*vertex)[3] = mesh->vertex;
//	float (*normal)[3] = mesh->normal;
//	int **link = mesh->vertex_link_v;
//	int *degree = mesh->degree_v;
//
//	float noise;
//	srand(100);
//	float size = mesh->averageOfEdgeLength();
//	int i;  for(i = 0; i<vertex_N; i++){
//		double r1, r2, r;
//		do{
//			r1 = 2.0*((double)rand())/RAND_MAX - 1.0;
//			r2 = 2.0*((double)rand())/RAND_MAX - 1.0;
//			r = r1*r1 + r2*r2;
//		}while(r >= 1.0 || r == 0.0);
//		double f = sqrt(-2.0*log(r)/r);
//		
//if(vertex[i][1] > 0)
//r1 = 0;
//		noise = (float)(size*rate*r1*f);
//		vertex[i][0] += normal[i][0]*noise;
//		vertex[i][1] += normal[i][1]*noise;
//		vertex[i][2] += normal[i][2]*noise;
//
//		if(++i == vertex_N)
//			return;
//		else{
//			
//if(vertex[i][1] > 0)
//r2 = 0;
//			noise = (float)(size*rate*r2*f);		
//			vertex[i][0] += normal[i][0]*noise;
//			vertex[i][1] += normal[i][1]*noise;
//			vertex[i][2] += normal[i][2]*noise;
//		}
//	}
//}

//void Smoother::checkGauusianSupportDikstra(int f, float sigma, float w)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int(*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//
//	int i;  for(i = 0; i<face_N; i++){
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0){
//				//dist[i][j] = 100000000;
//				continue;
//			}
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//			}
//		}
//	}
//
//	int *id = mesh->triangle_id;
//	if(mesh->triangle_id == NULL){
//		mesh->triangle_id = new int[face_N];
//		id = mesh->triangle_id;
//		for(i=0; i<face_N; i++)
//			id[i] = -1;
//	}
//
//	float *nf = normal[f];
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	bool *visit = new bool[face_N];
//	for(i=0; i<face_N; i++)
//		visit[i] = false;
//	geo[f] = 0;
//	int u = f;
//	while(true){
//		id[u] = (int)(255 - 255*geo[u]/(4*sigma));
//
//		int *nei = link[u];
//		int j;  for(i = 0; j<3; j++){
//			int v = nei[j];
//			if(v < 0 || id[v] >= 0)
//				continue;
//
//			double w1 = 0;//(1.0 - MeshData::DOT(nf,normal[u]));
//			double w2 = 0; //(1.0 - MeshData::DOT(nf,normal[v]));
//			double w3 = (1.0 - MeshData::DOT(normal[u],normal[v])); //3.0;
//			double g = (1.0+10*w*(w1+w2+w3))*dist[u][j];	
//			double M = geo[u] + g;
//			if(M > 4*sigma)
//				continue;
//			if(!visit[v]){
//				geo[v] = M;
//				//insert;
//				heap[++last_heap] = v;
//				index[v] = last_heap;
//				visit[v] = true;
//				//upheap;
//				upheap(geo, last_heap, last_heap, heap, index);
//			}
//			else if(M < geo[v]){
//				geo[v] = M;
//				//change;
//				if(index[v] != 1 && M < geo[heap[index[v]/2]])
//					//upheap;
//					upheap(geo, last_heap, index[v], heap, index);
//				else
//					//downheap;
//					downheap(geo, last_heap, index[v], heap, index);
//			}
//		}
//		if(last_heap == 0)
//			break;		
//
//		//delete;
//		u = heap[1];
//		heap[1] = heap[last_heap--];
//		index[heap[1]] = 1;
//		//down heap;
//		downheap(geo, last_heap, 1, heap, index);
//		if(last_heap == 0)
//			break;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//}

inline void Smoother::upheap(double *a, int N, int k, int *p, int *q)
{
	int v;
	v = p[k];
	while(k > 1 && a[p[k/2]] >= a[v]){
		p[k] = p[k/2]; q[p[k/2]] = k; k = k/2;
	}
	p[k] = v; q[v] = k;
}

//inline void Smoother::downheap(double *a, int N, int k, int *p, int *q)
//{
//	int j, v;
//	v = p[k];
//	while(k <= N/2){
//		j = k+k;
//		if(j < N && a[p[j]] > a[p[j+1]]) j++;
//		if(a[v] <= a[p[j]]) break;
//		p[k] = p[j]; q[p[j]] = k; k = j;
//	}
//	p[k] = v; q[v] = k;
//}

//void Smoother::computeOptimalGaussianDikstra(float s_min, float s_max, float step, float c)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//	int *visit = new int[face_N];
//	bool *decide = new bool[face_N];
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//	//float c = 0.001;
//	//float c = 0.00005;
//
//	//float total_e = 0;
//	//int n = 0;
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//				//n++;
//				//total_e += dist[i][j];
//			}
//		}
//		visit[i] = -1;
//		decide[i] = false;
//	}
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	
//	for(i=0; i<face_N; i++){
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//		}
//
//		geo[i] = 0;
//		int u = i;
//		visit[i] = i;
//		while(true){
//			decide[u] = true;
//
//			for(n=0; n<N; n++){
//				if(geo[u] > 4*si[n])
//					continue;
//				double w;
//				if(geo[u] <= 2.0*si[n])
//					w = area[u]*exp(-geo[u]*geo[u]/(2*si[n]*si[n])); 
//				else
//					w = area[u]*0.0625*pow(4.0 - geo[u]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[u][0];
//				ns[n][1] += w*normal[u][1];
//				ns[n][2] += w*normal[u][2];
//
//				total_w[n] += w;
//			}
//
//			int *nei = link[u];
//			int j;  for(i = 0; j<3; j++){
//				int v = nei[j];
//				if(v < 0 || (visit[v] == i && decide[v]))
//					continue;
//
//				double M = geo[u] + dist[u][j];
//				if(M > 4*s_max)
//					continue;
//				if(visit[v] != i){
//					geo[v] = M;
//					//insert;
//					heap[++last_heap] = v;
//					index[v] = last_heap;
//					visit[v] = i;
//					decide[v] = false;
//					//upheap;
//					upheap(geo, last_heap, last_heap, heap, index);
//				}
//				else if(M < geo[v]){
//					geo[v] = M;
//					//change;
//					if(index[v] != 1 && M < geo[heap[index[v]/2]])
//						//upheap;
//						upheap(geo, last_heap, index[v], heap, index);
//					else
//						//downheap;
//						downheap(geo, last_heap, index[v], heap, index);
//				}
//			}
//			if(last_heap == 0)
//				break;		
//
//			//delete;
//			u = heap[1];
//			heap[1] = heap[last_heap--];
//			index[heap[1]] = 1;
//			//down heap;
//			downheap(geo, last_heap, 1, heap, index);
//			if(last_heap == 0)
//				break;
//		}
//		
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			/*
//			double len = MeshData::LENGTH(ns[n]);
//			if((float)len == 0)
//				continue;
//			ns[n][0] /= len;
//			ns[n][1] /= len;
//			ns[n][2] /= len;
//			*/
//			
//			//double v = c/(si[n]*si[n]) + 2.0*(1.0-MeshData::DOT(ni, ns[n]));
//			double v = c/(si[n]*si[n]) + 1.0-MeshData::LENGTH(ns[n])/total_w[n];
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns[opt]);
//		
//		tmp_normal[i][0] = ns[opt][0]/len;
//		tmp_normal[i][1] = ns[opt][1]/len;
//		tmp_normal[i][2] = ns[opt][2]/len;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//	delete[] decide;
//	delete[] area;
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//
//	for(i=0; i<face_N; i++){
//		normal[i][0] = (float)tmp_normal[i][0];
//		normal[i][1] = (float)tmp_normal[i][1];
//		normal[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//}

//void Smoother::computeAnistoropicGaussianDikstra(float s_min, float s_max, float step, float c, float T)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//	int *visit = new int[face_N];
//	bool *decide = new bool[face_N];
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//	//float c = 0.001;
//	//float c = 0.00005;
//
//	//float total_e = 0;
//	//int n = 0;
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//				//n++;
//				//total_e += dist[i][j];
//			}
//		}
//		visit[i] = -1;
//		decide[i] = false;
//	}
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	
//	for(i=0; i<face_N; i++){
//		float *nf = normal[i];
//
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//		}
//		geo[i] = 0;
//		int u = i;
//		visit[i] = i;
//		while(true){
//			decide[u] = true;
//
//			for(n=0; n<N; n++){
//				if(geo[u] > 4*si[n])
//					continue;
//				double w;
//				if(geo[u] <= 2.0*si[n])
//					w = area[u]*exp(-geo[u]*geo[u]/(2*si[n]*si[n])); 
//				else
//					w = area[u]*0.0625*pow(4.0 - geo[u]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[u][0];
//				ns[n][1] += w*normal[u][1];
//				ns[n][2] += w*normal[u][2];
//
//				total_w[n] += w;
//			}
//
//			int *nei = link[u];
//			int j;  for(i = 0; j<3; j++){
//				int v = nei[j];
//				if(v < 0 || (visit[v] == i && decide[v]))
//					continue;
//
//				double w1 = (1.0 - MeshData::DOT(nf,normal[u]));
//				double w2 = (1.0 - MeshData::DOT(nf,normal[v]));
//				double w3 = (1.0 - MeshData::DOT(normal[u],normal[v]))/3.0;
//				double g = (1.0+T*(w1+w2+w3))*dist[u][j];	
//				double M = geo[u] + g;
//				if(M > 4*s_max)
//					continue;
//				if(visit[v] != i){
//					geo[v] = M;
//					//insert;
//					heap[++last_heap] = v;
//					index[v] = last_heap;
//					visit[v] = i;
//					decide[v] = false;
//					//upheap;
//					upheap(geo, last_heap, last_heap, heap, index);
//				}
//				else if(M < geo[v]){
//					geo[v] = M;
//					//change;
//					if(index[v] != 1 && M < geo[heap[index[v]/2]])
//						//upheap;
//						upheap(geo, last_heap, index[v], heap, index);
//					else
//						//downheap;
//						downheap(geo, last_heap, index[v], heap, index);
//				}
//			}
//			if(last_heap == 0)
//				break;		
//
//			//delete;
//			u = heap[1];
//			heap[1] = heap[last_heap--];
//			index[heap[1]] = 1;
//			//down heap;
//			downheap(geo, last_heap, 1, heap, index);
//			if(last_heap == 0)
//				break;
//		}
//		
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			/*
//			double len = MeshData::LENGTH(ns[n]);
//			if((float)len == 0)
//				continue;
//			ns[n][0] /= len;
//			ns[n][1] /= len;
//			ns[n][2] /= len;
//			*/
//			//double v = c/(si[n]*si[n]) + 2.0*(1.0-MeshData::DOT(ni, ns[n]));
//			double v = c/(si[n]*si[n]) + 1.0-MeshData::LENGTH(ns[n])/total_w[n];
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns[opt]);
//		
//		tmp_normal[i][0] = ns[opt][0]/len;
//		tmp_normal[i][1] = ns[opt][1]/len;
//		tmp_normal[i][2] = ns[opt][2]/len;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//	delete[] decide;
//	delete[] area;
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//
//	for(i=0; i<face_N; i++){
//		normal[i][0] = (float)tmp_normal[i][0];
//		normal[i][1] = (float)tmp_normal[i][1];
//		normal[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//}

//void Smoother::computeAnistoropicGaussianDikstraAngle(float s_min, float s_max, float step, float c, float T)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//	int *visit = new int[face_N];
//	bool *decide = new bool[face_N];
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//			}
//		}
//		visit[i] = -1;
//		decide[i] = false;
//	}
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	
//	for(i=0; i<face_N; i++){
//		float *nf = normal[i];
//
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//		}
//		geo[i] = 0;
//		int u = i;
//		visit[i] = i;
//		Node* Q = new Node;
//		while(true){
//			Q->append(-1, u);
//			decide[u] = true;
//
//			for(n=0; n<N; n++){
//				if(geo[u] > 4*si[n])
//					continue;
//				double w;
//				if(geo[u] <= 2.0*si[n])
//					w = area[u]*exp(-geo[u]*geo[u]/(si[n]*si[n])); 
//				else
//					w = area[u]*0.0625*pow(4.0 - geo[u]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[u][0];
//				ns[n][1] += w*normal[u][1];
//				ns[n][2] += w*normal[u][2];
//
//				total_w[n] += w;
//			}
//
//			int *nei = link[u];
//			int j;  for(i = 0; j<3; j++){
//				int v = nei[j];
//				if(v < 0 || (visit[v] == i && decide[v]))
//					continue;
//
//				double w1 = (1.0 - MeshData::DOT(nf,normal[u]));
//				double w2 = (1.0 - MeshData::DOT(nf,normal[v]));
//				double w3 = (1.0 - MeshData::DOT(normal[u],normal[v]))/3.0;
//				double g = (1.0+T*(w1+w2+w3))*dist[u][j];	
//				double M = geo[u] + g;
//				if(M > 4*s_max)
//					continue;
//				if(visit[v] != i){
//					geo[v] = M;
//					//insert;
//					heap[++last_heap] = v;
//					index[v] = last_heap;
//					visit[v] = i;
//					decide[v] = false;
//					//upheap;
//					upheap(geo, last_heap, last_heap, heap, index);
//				}
//				else if(M < geo[v]){
//					geo[v] = M;
//					//change;
//					if(index[v] != 1 && M < geo[heap[index[v]/2]])
//						//upheap;
//						upheap(geo, last_heap, index[v], heap, index);
//					else
//						//downheap;
//						downheap(geo, last_heap, index[v], heap, index);
//				}
//			}
//			if(last_heap == 0)
//				break;		
//
//			//delete;
//			u = heap[1];
//			heap[1] = heap[last_heap--];
//			index[heap[1]] = 1;
//			//down heap;
//			downheap(geo, last_heap, 1, heap, index);
//			if(last_heap == 0)
//				break;
//		}
//		
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			double len = MeshData::LENGTH(ns[n]);
//			if((float)len == 0)
//				continue;
//			ns[n][0] /= len;
//			ns[n][1] /= len;
//			ns[n][2] /= len;
//			
//			double angle = 0;
//			for(Node* current=Q; current->next!=NULL; current=current->next){
//				int c = current->f;
//				if(geo[c] > 4.0f*si[n])
//					continue;
//				double w;
//				if(geo[c] <= 2.0*si[n])
//					w = area[c]*exp(-geo[c]*geo[c]/(si[n]*si[n])); 
//				else
//					w = area[c]*0.0625*pow(4.0 - geo[c]/si[n], 4.0)/(E*E);
//				double dot = MeshData::DOT(ns[n], normal[c]);
//				angle += w*acos(dot);
//			}
//			double v = c/(si[n]*si[n]) + angle/(PI*total_w[n]);;
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns[opt]);
//		
//		tmp_normal[i][0] = ns[opt][0]; //len;
//		tmp_normal[i][1] = ns[opt][1]; //len;
//		tmp_normal[i][2] = ns[opt][2]; //len;
//
//		delete Q;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//	delete[] decide;
//	delete[] area;
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//
//	for(i=0; i<face_N; i++){
//		normal[i][0] = (float)tmp_normal[i][0];
//		normal[i][1] = (float)tmp_normal[i][1];
//		normal[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//}

//void Smoother::checkGaussianSupportDikstra2(int f, float sigma, float w)
//{
//	int face_N = mesh->face_N;
//	float (*normal_o)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int(*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//
//	mesh->normal_f = NULL;
//	mesh->computeFaceNormal();
//	//int i;  for(i = 0; i<5; i++)
//		//smoothNormal(0, 0, 0);
//	float (*normal)[3] = mesh->normal_f;
//
//	int i;  for(i = 0; i<face_N; i++){
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0){
//				//dist[i][j] = 100000000;
//				continue;
//			}
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//
//				double PM = MeshData::DIST(center[i], p);
//				double MQ = MeshData::DIST(p, center[ad]);
//				//float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				double dot = MeshData::DOT(normal[i], normal[ad]);
//				dist[i][j] = (float)(PM + MQ + 
//					w*(PM*PM + MQ*MQ -2.0*dot*PM*MQ)/pow(PM+MQ, 3.0));
//
//				/*
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				//float dn = w*2.0*(1.0 - MeshData::DOT(normal[i], normal[ad]));
//				double dot = MeshData::DOT(normal[i], normal[ad]);
//				if(dot > 1.0)
//					dot = 1;
//				else if(dot < -1.0)
//					dot = -1;
//				double a = acos(dot);
//				dist[i][j] = (1.0f + (float)(w*a*a))*d;
//				*/
//			}
//		}
//	}
//
//	int *id = mesh->triangle_id;
//	if(mesh->triangle_id == NULL){
//		mesh->triangle_id = new int[face_N];
//		id = mesh->triangle_id;
//		for(i=0; i<face_N; i++)
//			id[i] = -1;
//	}
//
//	float *nf = normal[f];
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	bool *visit = new bool[face_N];
//	for(i=0; i<face_N; i++)
//		visit[i] = false;
//	geo[f] = 0;
//	int u = f;
//	while(true){
//		id[u] = (int)(255 - 255*geo[u]/(4*sigma));
//
//		int *nei = link[u];
//		int j;  for(i = 0; j<3; j++){
//			int v = nei[j];
//			if(v < 0 || id[v] >= 0)
//				continue;
//
//			double M = geo[u] + dist[u][j];
//			if(M > 4*sigma)
//				continue;
//			if(!visit[v]){
//				geo[v] = M;
//				//insert;
//				heap[++last_heap] = v;
//				index[v] = last_heap;
//				visit[v] = true;
//				//upheap;
//				upheap(geo, last_heap, last_heap, heap, index);
//			}
//			else if(M < geo[v]){
//				geo[v] = M;
//				//change;
//				if(index[v] != 1 && M < geo[heap[index[v]/2]])
//					//upheap;
//					upheap(geo, last_heap, index[v], heap, index);
//				else
//					//downheap;
//					downheap(geo, last_heap, index[v], heap, index);
//			}
//		}
//		if(last_heap == 0)
//			break;		
//
//		//delete;
//		u = heap[1];
//		heap[1] = heap[last_heap--];
//		index[heap[1]] = 1;
//		//down heap;
//		downheap(geo, last_heap, 1, heap, index);
//		if(last_heap == 0)
//			break;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//
//	delete[] normal;
//	mesh->normal_f = normal_o;
//}

//void Smoother::computeOptimalGaussianDikstra2(float s_min, float s_max, float step, float c)
//{
//	int face_N = mesh->face_N;
//	float (*normal)[3] = mesh->normal_f;
//	mesh->normal_f = NULL;
//	mesh->computeFaceNormal();
//	float (*normal_o)[3] = mesh->normal_f;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//	int *visit = new int[face_N];
//	bool *decide = new bool[face_N];
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//	//float c = 0.001;
//	//float c = 0.00005;
//
//	//float total_e = 0;
//	//int n = 0;
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				dist[i][j] = d;
//				//n++;
//				//total_e += dist[i][j];
//			}
//		}
//		visit[i] = -1;
//		decide[i] = false;
//	}
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	double (*ns_o)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	
//	for(i=0; i<face_N; i++){
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//			ns_o[n][0] = 0;
//			ns_o[n][1] = 0;
//			ns_o[n][2] = 0;
//		}
//
//		geo[i] = 0;
//		int u = i;
//		visit[i] = i;
//		while(true){
//			decide[u] = true;
//
//			for(n=0; n<N; n++){
//				if(geo[u] > 4*si[n])
//					continue;
//				double w;
//				if(geo[u] <= 2.0*si[n])
//					w = area[u]*exp(-geo[u]*geo[u]/(2*si[n]*si[n])); 
//				else
//					w = area[u]*0.0625*pow(4.0 - geo[u]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[u][0];
//				ns[n][1] += w*normal[u][1];
//				ns[n][2] += w*normal[u][2];
//
//				ns_o[n][0] += w*normal_o[u][0];
//				ns_o[n][1] += w*normal_o[u][1];
//				ns_o[n][2] += w*normal_o[u][2];
//
//				total_w[n] += w;
//			}
//
//			int *nei = link[u];
//			int j;  for(i = 0; j<3; j++){
//				int v = nei[j];
//				if(v < 0 || (visit[v] == i && decide[v]))
//					continue;
//
//				double M = geo[u] + dist[u][j];
//				if(M > 4*s_max)
//					continue;
//				if(visit[v] != i){
//					geo[v] = M;
//					//insert;
//					heap[++last_heap] = v;
//					index[v] = last_heap;
//					visit[v] = i;
//					decide[v] = false;
//					//upheap;
//					upheap(geo, last_heap, last_heap, heap, index);
//				}
//				else if(M < geo[v]){
//					geo[v] = M;
//					//change;
//					if(index[v] != 1 && M < geo[heap[index[v]/2]])
//						//upheap;
//						upheap(geo, last_heap, index[v], heap, index);
//					else
//						//downheap;
//						downheap(geo, last_heap, index[v], heap, index);
//				}
//			}
//			if(last_heap == 0)
//				break;		
//
//			//delete;
//			u = heap[1];
//			heap[1] = heap[last_heap--];
//			index[heap[1]] = 1;
//			//down heap;
//			downheap(geo, last_heap, 1, heap, index);
//			if(last_heap == 0)
//				break;
//		}
//		
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			double v = c/(si[n]*si[n]) + 1.0-MeshData::LENGTH(ns[n])/total_w[n];
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns_o[opt]);
//		
//		tmp_normal[i][0] = ns_o[opt][0]/len;
//		tmp_normal[i][1] = ns_o[opt][1]/len;
//		tmp_normal[i][2] = ns_o[opt][2]/len;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//	delete[] decide;
//	delete[] area;
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//	delete[] normal;
//	delete[] ns_o;
//
//	for(i=0; i<face_N; i++){
//		normal_o[i][0] = (float)tmp_normal[i][0];
//		normal_o[i][1] = (float)tmp_normal[i][1];
//		normal_o[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//}

//void Smoother::computeAnistoropicGaussianDikstra2(float s_min, float s_max, float step, float c, float T)
//{
//	int face_N = mesh->face_N;
//	int (*link)[3] = mesh->face_link_E;
//	double (*tmp_normal)[3] = new double[face_N][3];
//	double *area = new double[face_N];
//	double *geo = new double[face_N];
//	float (*dist)[3] = new float[face_N][3];
//	mesh->computeCenter();
//	float (*center)[3] = mesh->center;
//	int (*face)[3] = mesh->face;
//	float (*vertex)[3] = mesh->vertex;
//	int *visit = new int[face_N];
//	bool *decide = new bool[face_N];
//
//	if(sigma != NULL)
//		delete[] sigma;
//	sigma = new float[face_N];
//	//float c = 0.001;
//	//float c = 0.00005;
//
//	//float total_e = 0;
//	//int n = 0;
//
//	//int i;  for(i = 0; i<5; i++)
//		//smoothNormal(0,0,0);
//	float (*normal)[3] = mesh->normal_f;
//
//	int i;  for(i = 0; i<face_N; i++){
//		area[i] = mesh->faceArea(i);
//		int j;  for(i = 0; j<3; j++){
//			int ad = link[i][j];
//			if(ad < 0)
//				continue;
//			if(ad < i){
//				int k;
//				if(link[ad][0] == i)
//					k = 0;
//				else if(link[ad][1] == i)
//					k = 1;
//				else
//					k = 2;
//				dist[i][j] = dist[ad][k];
//			}
//			else{
//				float *p1 = vertex[face[i][j]];
//				float *p2 = vertex[face[i][(j+1)%3]];
//				float p[3];
//				p[0] = 0.5f*(p1[0] + p2[0]);
//				p[1] = 0.5f*(p1[1] + p2[1]);
//				p[2] = 0.5f*(p1[2] + p2[2]);
//				//float d1 = (float)MeshData::DIST(center[i], center[ad]);
//				double PM = MeshData::DIST(center[i], p);
//				double MQ = MeshData::DIST(p, center[ad]);
//				//float d =  (float)(MeshData::DIST(center[i], p) + MeshData::DIST(p, center[ad]));
//				double dot = MeshData::DOT(normal[i], normal[ad]);
//				dist[i][j] = (float)(PM + MQ + 
//					T*(PM*PM + MQ*MQ -2.0*dot*PM*MQ)/pow(PM+MQ, 3.0));
//				/*
//				if(dot > 1.0)
//					dot = 1;
//				else if(dot < -1.0)
//					dot = -1;
//				double a = acos(dot);
//				*/
//				//dist[i][j] = (1.0f + (float)(T*a*a))*d;
//				//dist[i][j] = d;
//				//n++;
//				//total_e += dist[i][j];
//			}
//		}
//		visit[i] = -1;
//		decide[i] = false;
//	}
//
//	mesh->computeFaceNormal();
//	normal = mesh->normal_f;
//
//	int N = (int)((s_max-s_min)/step) + 1;
//	double (*ns)[3] = new double[N][3];
//	float *si = new float[N];
//	double *total_w = new double[N];
//	for(i=0; i<N; i++)
//		si[i] = s_min + step*i;
//
//	int* heap = new int[face_N];
//	int* index = new int[face_N];
//	int last_heap = 0;
//	
//	for(i=0; i<face_N; i++){
//		float *nf = normal[i];
//
//		for(int n=0; n<N; n++){
//			ns[n][0] = 0;
//			ns[n][1] = 0;
//			ns[n][2] = 0;
//			total_w[n] = 0;
//		}
//		geo[i] = 0;
//		int u = i;
//		visit[i] = i;
//		while(true){
//			decide[u] = true;
//
//			for(n=0; n<N; n++){
//				if(geo[u] > 4*si[n])
//					continue;
//				double w;
//				if(geo[u] <= 2.0*si[n])
//					w = area[u]*exp(-geo[u]*geo[u]/(2*si[n]*si[n])); 
//				else
//					w = area[u]*0.0625*pow(4.0 - geo[u]/si[n], 4.0)/(E*E);
//
//				ns[n][0] += w*normal[u][0];
//				ns[n][1] += w*normal[u][1];
//				ns[n][2] += w*normal[u][2];
//
//				total_w[n] += w;
//			}
//
//			int *nei = link[u];
//			int j;  for(i = 0; j<3; j++){
//				int v = nei[j];
//				if(v < 0 || (visit[v] == i && decide[v]))
//					continue;
//
//				double M = geo[u] + dist[u][j];
//				if(M > 4*s_max)
//					continue;
//				if(visit[v] != i){
//					geo[v] = M;
//					//insert;
//					heap[++last_heap] = v;
//					index[v] = last_heap;
//					visit[v] = i;
//					decide[v] = false;
//					//upheap;
//					upheap(geo, last_heap, last_heap, heap, index);
//				}
//				else if(M < geo[v]){
//					geo[v] = M;
//					//change;
//					if(index[v] != 1 && M < geo[heap[index[v]/2]])
//						//upheap;
//						upheap(geo, last_heap, index[v], heap, index);
//					else
//						//downheap;
//						downheap(geo, last_heap, index[v], heap, index);
//				}
//			}
//			if(last_heap == 0)
//				break;		
//
//			//delete;
//			u = heap[1];
//			heap[1] = heap[last_heap--];
//			index[heap[1]] = 1;
//			//down heap;
//			downheap(geo, last_heap, 1, heap, index);
//			if(last_heap == 0)
//				break;
//		}
//		
//		float *ni = normal[i];
//		double min = 100000000;
//		int opt;
//		for(n=0; n<N; n++){
//			/*
//			double len = MeshData::LENGTH(ns[n]);
//			if((float)len == 0)
//				continue;
//			ns[n][0] /= len;
//			ns[n][1] /= len;
//			ns[n][2] /= len;
//			*/
//			//double v = c/(si[n]*si[n]) + 2.0*(1.0-MeshData::DOT(ni, ns[n]));
//			double v = c/(si[n]*si[n]) + 1.0-MeshData::LENGTH(ns[n])/total_w[n];
//			if(min > v){
//				opt = n;
//				min = v;
//			}
//		}
//		sigma[i] = (si[opt] - s_min)/(s_max - s_min);
//		double len = MeshData::LENGTH(ns[opt]);
//		
//		tmp_normal[i][0] = ns[opt][0]/len;
//		tmp_normal[i][1] = ns[opt][1]/len;
//		tmp_normal[i][2] = ns[opt][2]/len;
//	}
//	delete[] index;
//	delete[] heap;
//	delete[] dist;
//	delete[] geo;
//	delete[] visit;
//	delete[] decide;
//	delete[] area;
//	delete[] ns;
//	delete[] si;
//	delete[] total_w;
//
//	for(i=0; i<face_N; i++){
//		normal[i][0] = (float)tmp_normal[i][0];
//		normal[i][1] = (float)tmp_normal[i][1];
//		normal[i][2] = (float)tmp_normal[i][2];
//	}
//
//	delete[] tmp_normal;
//}
