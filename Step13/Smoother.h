/************************************************************************
Yutaka Ohtake 

Smoother.h
Smoothing methods

Copyright (c) 1999-2001 
The University of Aizu. All Rights Reserved.
************************************************************************/

#pragma once


#include "MeshData.h"
#include "Node.h"
//#include "PBCG.h"

#define EPSILON 0.00001

class Smoother  
{
public:
	MeshData* mesh;
	float (*T)[3];
	float *sigma;

public:
	void upheap(double *a, int N, int k, int *p, int *q);
//	void downheap(double *a, int N, int k, int *p, int *q);
	void upheap(float *a, int N, int k, int *p, int *q);
//	void downheap(float *a, int N, int k, int *p, int *q);
	void smoothRidgeDir(int iter, float T);
//	void moveMinPosition();
	void smoothNormalV(int iter, float T, float rate);
	void smoothNormalV(float T, float rate, int f_type, int d_type);
//	void curvatureAlongMedian(int iter, float dt);
//	void adaptiveSmoothing(int iter, float T, float rate, int int_type, int inner_iter, float int_step, int d_type, int ave_iter, int func_type, int L_deg);
//	void moveNormal(float dt);
//	void smoothNormal(float T, int f_type, int d_type);
//	void smoothTagEdges(int iter, float dt);
//	void smoothKmaxKmin(int iter, float T);
//	void smoothTmaxTmin(int iter, float T);
//	void LaplacianWithRidge(int iter, float dt);
//	void meanCurvatureFlow12Ring(int iter, float dt, float e);
//	void computeNRingMeanCurvature(int n);
//	void GaussianNoise(float size);
	void setMeshData(MeshData* mesh);
//	void meanCurvatureFlow(int iter, float dt);
//	void LaplacianFlow(int iter, float dt);
	Smoother();
	virtual ~Smoother();

	static inline BOOL INVERSE(double B[6], double A[6]){
		double d = DET(A);
		if(fabs(d) < EPSILON)
			return false;
		B[0] = (A[3]*A[5] - A[4]*A[4])/d;
		B[1] = (A[2]*A[4] - A[1]*A[5])/d;
		B[2] = (A[1]*A[4] - A[2]*A[3])/d;
		B[3] = (A[0]*A[5] - A[2]*A[2])/d;
		B[4] = (A[1]*A[2] - A[0]*A[4])/d;
		B[5] = (A[0]*A[3] - A[1]*A[1])/d;
		return true;
	}

	static inline double DET(double A[6]){
		return A[0]*A[3]*A[5] + 2.0*A[1]*A[4]*A[2] 
			-A[2]*A[2]*A[3] - A[1]*A[1]*A[5] - A[4]*A[4]*A[0];
	}

	static inline void MATRIX(double A[6], double n[3]){
		A[0] = n[0]*n[0];
		A[1] = n[0]*n[1];
		A[2] = n[0]*n[2];
		A[3] = n[1]*n[1];
		A[4] = n[1]*n[2];
		A[5] = n[2]*n[2];
	}

	static inline void MAT_TIMES(double A[6], double k){
		A[0] *= k;
		A[1] *= k;
		A[2] *= k;
		A[3] *= k;
		A[4] *= k;
		A[5] *= k;
	}

	static inline void MAT_SUM(double B[6], double A[6]){
		B[0] += A[0];
		B[1] += A[1];
		B[2] += A[2];
		B[3] += A[3];
		B[4] += A[4];
		B[5] += A[5];
	}

	static inline void MAT_BY_VEC(double v[3], double A[6], double b[3]){
		v[0] = A[0]*b[0] + A[1]*b[1] + A[2]*b[2];
		v[1] = A[1]*b[0] + A[3]*b[1] + A[4]*b[2];
		v[2] = A[2]*b[0] + A[4]*b[1] + A[5]*b[2];
	}

	public:
//		void computeAnistoropicGaussianDikstra2(float s_min, float s_max, float step, float c, float T);
//		void computeOptimalGaussianDikstra2(float s_min, float s_max, float step, float c);
//		void checkGaussianSupportDikstra2(int f, float sigma, float w);
//		void computeAnistoropicGaussianDikstraAngle(float s_min, float s_max, float step, float c, float T);
//		void computeAnistoropicGaussianDikstra(float s_min, float s_max, float step, float c, float T);
//		void computeOptimalGaussianDikstra(float s_min, float s_max, float step, float c);
//		void checkGauusianSupportDikstra(int f, float sigma, float w);
//		void addEdgeLengthNoise(float rate);
//		void computeAnistoropicGaussian(float s_min, float s_max, float step, float c, float T);
//		void adaptiveGaussianSmoothing(float c, BOOL is_anistoropic, float w);
//		float computeAliasingMeasure(int i1, int i2);
//		void antialiasingEdgeFlip();
//		bool flipEdge(int i1, int i2, int &i3, int &i4);
//		void computeOptimalGaussian(float s_min, float s_max, float step, float c);
//		float computeSigma2(int f, int size);
//		float computeSigma1(int f, int size);
//		void checkGaussianSupport(int f, float sigma, float T);
//		void decideSigma(float *sigma, float c, int size);
//		void smoothNormalGaussian(int size, float *sigma);
//		void minimizeNormalErr();
//		void IntegrateNormalImplicit(float dt);
//		void MeanCurvatureFlowImplicit(float dt);
//		void setupPBCG(PBCG* pbcg);
//		void LaplacianFlowImplicit(float dt);
//		void adaptiveSmoothingM(int iter, BOOL isWeighted, int int_type, int inner_iter, float int_step, int d_type, int ave_iter, int L_deg);
		void smoothNormalMedian(BOOL isWeighted, int d_type);
//		void badTriangles(BOOL *isBad);
		void smoothNormalMori(int iter);
//		void Bilaplacian(int iter, float dt);
//		void connectTag(int iter, float dt, float angle);
//		void smoothRidgeEdges(int iter, float dt);
//		void moveNormalAreaD(float dt, int power);
		void smoothNormalV2(float T, float rate, int f_type);
//		void moveNormalAreaD(float dt);
		void TaubinMethod(int iter, float dt, float low_pass, int type);
//		void LaplacianLocalControl(int iter, float dt);
//		void moveNormalTaubin(float C);
//		void BilaplacianWithRidge2(int iter, float dt);
//		void BilaplacianWithRidge(int iter, float dt);
//		void moveNormalD(float dt);
//		void smoothTmaxTminSIG2();
//		void smoothTmaxTminSIG1(int iter, float T);
//		void smoothNormalL();
//		void smoothTmaxTmin4(int iter, float T);
		void smoothKmaxKmin2(int iter, float T);
//		void smoothTmaxTmin3(int iter, float T);
		void smoothRidgeDir3(int iter, float T);
		void smoothRidgeDir2(int iter, float T);
//		void smoothTmaxTmin2(int iter, float T);

	static inline void MAT_VEC(double y[3], double A[3][3], double x[3]){
		y[0] = A[0][0]*x[0] + A[0][1]*x[1] + A[0][2]*x[2];
		y[1] = A[1][0]*x[0] + A[1][1]*x[1] + A[1][2]*x[2];
		y[2] = A[2][0]*x[0] + A[2][1]*x[1] + A[2][2]*x[2];
	}

	static inline void MAT_TIMES(double C[3][3], double A[3][3], double B[3][3]){
		C[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
		C[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
		C[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];

		C[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
		C[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
		C[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];

		C[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
		C[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
		C[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
	}

	static inline void GENERATE_MAT(double R[3][3], double angle, double v[3]){
		double l = sqrt(v[0]*v[0] + v[1]*v[1]);
		double a;
		if(l == 0)
			a = 0;
		else
			a = acos(v[1]/l);
		if(v[0] < 0)
			a = -a;
		double b = acos(v[2]/MeshData::LENGTH(v));

		double Rz[3][3], Rx[3][3], Rt[3][3]; 
		Rz[0][0] = Rz[1][1] = (float)cos(a);
		Rz[0][1] = sin(a);
		Rz[1][0] = -Rz[0][1];
		Rz[2][0] = Rz[2][1] = Rz[0][2] = Rz[1][2] = 0;
		Rz[2][2] = 1;

		Rt[0][0] = Rt[1][1] = cos(angle);
		Rt[0][1] = -sin(angle);
		Rt[1][0] = -Rt[0][1];
		Rt[2][0] = Rt[2][1] = Rt[0][2] = Rt[1][2] = 0;
		Rt[2][2] = 1;

		Rx[0][0] = 1;
		Rx[0][1] = Rx[0][2] = Rx[1][0] = Rx[2][0] = 0;
		Rx[1][1] = Rx[2][2] = cos(b);
		Rx[1][2] = sin(b);
		Rx[2][1] = -Rx[1][2];

		double tmp[3][3];
		MAT_TIMES(tmp, Rz, Rx);
		MAT_TIMES(R, tmp, Rt);
		Rz[1][0] *= -1;
		Rz[0][1] *= -1;
		Rx[1][2] *= -1;
		Rx[2][1] *= -1;
		MAT_TIMES(tmp, R, Rx);
		MAT_TIMES(R, tmp, Rz);

		MAT_TIMES(tmp, Rx, Rz);
		double x[3];
		MAT_VEC(x, tmp, v);
	}

	static inline void GENERATE_MAT(double R[3][3], double angle, float v1[3]){
		double v[3];
		v[0] = v1[0];
		v[1] = v1[1];
		v[2] = v1[2];

		double l = sqrt(v[0]*v[0] + v[1]*v[1]);
		double a;
		if(l == 0)
			a = 0;
		else
			a = acos(v[1]/l);
		if(v[0] < 0)
			a = -a;
		double b = acos(v[2]/MeshData::LENGTH(v));

		double Rz[3][3], Rx[3][3], Rt[3][3]; 
		Rz[0][0] = Rz[1][1] = (float)cos(a);
		Rz[0][1] = sin(a);
		Rz[1][0] = -Rz[0][1];
		Rz[2][0] = Rz[2][1] = Rz[0][2] = Rz[1][2] = 0;
		Rz[2][2] = 1;

		Rt[0][0] = Rt[1][1] = cos(angle);
		Rt[0][1] = -sin(angle);
		Rt[1][0] = -Rt[0][1];
		Rt[2][0] = Rt[2][1] = Rt[0][2] = Rt[1][2] = 0;
		Rt[2][2] = 1;

		Rx[0][0] = 1;
		Rx[0][1] = Rx[0][2] = Rx[1][0] = Rx[2][0] = 0;
		Rx[1][1] = Rx[2][2] = cos(b);
		Rx[1][2] = sin(b);
		Rx[2][1] = -Rx[1][2];

		double tmp[3][3];
		MAT_TIMES(tmp, Rz, Rx);
		MAT_TIMES(R, tmp, Rt);
		Rz[1][0] *= -1;
		Rz[0][1] *= -1;
		Rx[1][2] *= -1;
		Rx[2][1] *= -1;
		MAT_TIMES(tmp, R, Rx);
		MAT_TIMES(R, tmp, Rz);

		MAT_TIMES(tmp, Rx, Rz);
		double x[3];
		MAT_VEC(x, tmp, v);
	}

};

