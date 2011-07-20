/*
Copyright (c) 2011, Movania Muhammad Mobeen
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list
of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

//A simple cloth using implicit Euler integration based on the SIGGRAPH course notes
//"Realtime Physics" http://www.matthiasmueller.info/realtimephysics/coursenotes.pdf
//using GLUT,GLEW and GLM libraries. This code is intended for beginners
//so that they may understand what is required to implement implicit integration
//based cloth simulation.
//
//This code is under BSD license. If you make some improvements,
//or are using this in your research, do let me know and I would appreciate
//if you acknowledge this in your code.
//
//Controls:
//left click on any empty region to rotate, middle click to zoom
//left click and drag any point to drag it.
//
//Author: Movania Muhammad Mobeen
//        School of Computer Engineering,
//        Nanyang Technological University,
//        Singapore.
//Email : mova0002@e.ntu.edu.sg

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glut.h>
#include <vector>
#include <glm/glm.hpp>

#pragma comment(lib, "glew32.lib")

using namespace std;
const int width = 1024, height = 1024;

int numX = 20, numY=20;
const size_t total_points = (numX+1)*(numY+1);
int size = 4;
float hsize = size/2.0f;

char info[MAX_PATH]={0};

int selected_index = -1;

struct Spring {
	int p1, p2;
	float rest_length;
	float Ks, Kd;
	int type;
};

const float EPS = 0.001f;
const float EPS2 = EPS*EPS;
const int i_max = 10;

void StepPhysics(float dt);

void SolveConjugateGradient2(glm::mat2 A, glm::vec2& x, glm::vec2 b) {
	float i =0;
	glm::vec2 r = b - A*x;
	glm::vec2 d = r;
	glm::vec2 q = glm::vec2(0);
	float alpha_new = 0;
	float alpha = 0;
	float beta  = 0;
	float delta_old = 0;
	float delta_new = glm::dot(r,r);
	float delta0    = delta_new;
	while(i<i_max && delta_new> EPS2) {
		q = A*d;
		alpha = delta_new/glm::dot(d,q);
		x = x + alpha*d;
		r = r - alpha*q;
		delta_old = delta_new;
		delta_new = glm::dot(r,r);
		beta = delta_new/delta_old;
		d = r + beta*d;
		i++;
	}
}



template<class T>
class LargeVector {
private:
	vector<T> v;

public:

	LargeVector() {

	}
	LargeVector(const LargeVector& other) {
		v.resize(other.v.size());
		memcpy(&v[0], &(other.v[0]), sizeof(other.v[0])*other.v.size());
	}
	void resize(const int size) {
		v.resize(size);
	}
	void clear(bool isIdentity=false) {
		memset(&v[0], 0, sizeof(T)*v.size());
		if(isIdentity) {
			for(size_t i=0;i<v.size();i++) {
				v[i] = T(1);
			}
		}
	}
	size_t size() {
		return v.size();
	}


	T& operator[](int index) {
		return v[index];
	}
	friend LargeVector<glm::vec3> operator*(const LargeVector<glm::mat3> other, const LargeVector<glm::vec3> f );
	friend LargeVector<glm::vec3> operator*(const float f, const LargeVector<glm::vec3> other);
	friend LargeVector<glm::vec3> operator-(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb );
	//friend LargeVector<T> operator+(const LargeVector<T> Va, const LargeVector<T> Vb );
	friend LargeVector<glm::vec3> operator+(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb );

	friend LargeVector<glm::mat3> operator*(const float f, const LargeVector<glm::mat3> other);
	friend LargeVector<glm::mat3> operator-(const LargeVector<glm::mat3> Va, const LargeVector<glm::mat3> Vb );
	//friend LargeVector<glm::mat3> operator+(const LargeVector<glm::mat3> Va, const LargeVector<glm::mat3> Vb );


	friend LargeVector<glm::vec3> operator/(const float f, const LargeVector<glm::vec3> v );
	friend float dot(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb );
};

LargeVector<glm::vec3> operator*(const LargeVector<glm::mat3> other, const LargeVector<glm::vec3> v ) {
	LargeVector<glm::vec3> tmp(v);
	for(size_t i=0;i<v.v.size();i++) {
		tmp.v[i] = other.v[i] * v.v[i];
	}
	return tmp;
}

LargeVector<glm::vec3> operator*(const float f, const LargeVector<glm::vec3> other) {
	LargeVector<glm::vec3> tmp(other);
	for(size_t i=0;i<other.v.size();i++) {
		tmp.v[i] = other.v[i]*f;
	}
	return tmp;
}
LargeVector<glm::mat3> operator*(const float f, const LargeVector<glm::mat3> other) {
	LargeVector<glm::mat3> tmp(other);
	for(size_t i=0;i<other.v.size();i++) {
		tmp.v[i] = other.v[i]*f;
	}
	return tmp;
}
LargeVector<glm::vec3> operator-(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb ) {
	LargeVector<glm::vec3> tmp(Va);
	for(size_t i=0;i<Va.v.size();i++) {
		tmp.v[i] = Va.v[i] - Vb.v[i];
	}
	return tmp;
}
LargeVector<glm::mat3> operator-(const LargeVector<glm::mat3> Va, const LargeVector<glm::mat3> Vb ) {
	LargeVector<glm::mat3> tmp(Va);
	for(size_t i=0;i<Va.v.size();i++) {
		tmp.v[i] = Va.v[i] - Vb.v[i];
	}
	return tmp;
}

LargeVector<glm::vec3> operator+(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb ) {
	LargeVector<glm::vec3> tmp(Va);
	for(size_t i=0;i<Va.v.size();i++) {
		tmp.v[i] = Va.v[i] + Vb.v[i];
	}
	return tmp;
}

LargeVector<glm::vec3> operator/(const float f, const LargeVector<glm::vec3> v ) {
	LargeVector<glm::vec3> tmp(v);
	for(size_t i=0;i<v.v.size();i++) {
		tmp.v[i] = v.v[i] / f;
	}
	return tmp;
}


float dot(const LargeVector<glm::vec3> Va, const LargeVector<glm::vec3> Vb ) {
	float sum = 0;
	for(size_t i=0;i<Va.v.size();i++) {
		sum += glm::dot(Va.v[i],Vb.v[i]);
	}
	return sum;
}

void SolveConjugateGradient(LargeVector<glm::mat3> A, LargeVector<glm::vec3>& x, LargeVector<glm::vec3> b) {
	float i =0;
	LargeVector<glm::vec3> r = b - A*x;
	LargeVector<glm::vec3> d = r;
	LargeVector<glm::vec3> q;
	float alpha_new = 0;
	float alpha = 0;
	float beta  = 0;
	float delta_old = 0;
	float delta_new = dot(r,r);
	float delta0    = delta_new;
	while(i<i_max && delta_new> EPS2) {
		q = A*d;
		alpha = delta_new/dot(d,q);
		x = x + alpha*d;
		r = r - alpha*q;
		delta_old = delta_new;
		delta_new = dot(r,r);
		beta = delta_new/delta_old;
		d = r + beta*d;
		i++;
	}
}


void SolveConjugateGradient(glm::mat3 A, glm::vec3& x, glm::vec3 b) {
	float i =0;
	glm::vec3 r = b - A*x;
	glm::vec3 d = r;
	glm::vec3 q = glm::vec3(0);
	float alpha_new = 0;
	float alpha = 0;
	float beta  = 0;
	float delta_old = 0;
	float delta_new = glm::dot(r,r);
	float delta0    = delta_new;
	while(i<i_max && delta_new> EPS2) {
		q = A*d;
		alpha = delta_new/glm::dot(d,q);
		x = x + alpha*d;
		r = r - alpha*q;
		delta_old = delta_new;
		delta_new = glm::dot(r,r);
		beta = delta_new/delta_old;
		d = r + beta*d;
		i++;
	}
}




vector<GLushort> indices;
vector<Spring> springs;

LargeVector<glm::vec3> X;
LargeVector<glm::vec3> V;
LargeVector<glm::vec3> F;

LargeVector<glm::mat3> df_dx; //  df/dp
LargeVector<glm::vec3> dc_dp; //  df/dp

vector<glm::vec3> deltaP2;
LargeVector<glm::vec3> V_new;
LargeVector<glm::mat3> M; //the mass matrix
glm::mat3 I=glm::mat3(1);//identity matrix

vector<float> C; //for implicit integration
vector<float> C_Dot; //for implicit integration


int oldX=0, oldY=0;
float rX=15, rY=0;
int state =1 ;
float dist=-23;
const int GRID_SIZE=10;

const int STRUCTURAL_SPRING = 0;
const int SHEAR_SPRING = 1;
const int BEND_SPRING = 2;
int spring_count=0;


const float DEFAULT_DAMPING =  -0.0125f;
float	KsStruct = 0.75f,KdStruct = -0.25f;
float	KsShear = 0.75f,KdShear = -0.25f;
float	KsBend = 0.95f,KdBend = -0.25f;
glm::vec3 gravity=glm::vec3(0.0f,-0.00981f,0.0f);
float mass = 0.5f;

float timeStep =  1/60.0f;
float currentTime = 0;
double accumulator = timeStep;

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
double frameTimeQP=0;
float frameTime =0 ;


GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;
float startTime =0, fps=0;
int totalFrames=0;



void AddSpring(int a, int b, float ks, float kd, int type) {
	Spring spring;
	spring.p1=a;
	spring.p2=b;
	spring.Ks=ks;
	spring.Kd=kd;
	spring.type = type;
	glm::vec3 deltaP = X[a]-X[b];
	spring.rest_length = sqrt(glm::dot(deltaP, deltaP));
	springs.push_back(spring);
}
void OnMouseDown(int button, int s, int x, int y)
{
	if (s == GLUT_DOWN)
	{
		oldX = x;
		oldY = y;
		int window_y = (height - y);
		float norm_y = float(window_y)/float(height/2.0);
		int window_x = x ;
		float norm_x = float(window_x)/float(width/2.0);

		float winZ=0;
		glReadPixels( x, height-y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
		if(winZ==1)
			winZ=0;
		double objX=0, objY=0, objZ=0;
		gluUnProject(window_x,window_y, winZ,  MV,  P, viewport, &objX, &objY, &objZ);
		glm::vec3 pt(objX,objY, objZ);
		size_t i=0;
		for(i=0;i<total_points;i++) {
			if( glm::distance(X[i],pt)<0.1) {
				selected_index = i;
				printf("Intersected at %d\n",i);
				break;
			}
		}
	}

	if(button == GLUT_MIDDLE_BUTTON)
		state = 0;
	else
		state = 1;

	if(s==GLUT_UP) {
		selected_index= -1;
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}

void OnMouseMove(int x, int y)
{
	if(selected_index == -1) {
		if (state == 0)
			dist *= (1 + (y - oldY)/60.0f);
		else
		{
			rY += (x - oldX)/5.0f;
			rX += (y - oldY)/5.0f;
		}
	} else {
		float delta = 1500/abs(dist);
		float valX = (x - oldX)/delta;
		float valY = (oldY - y)/delta;
		if(abs(valX)>abs(valY))
			glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
		else
			glutSetCursor(GLUT_CURSOR_UP_DOWN);

		V[selected_index] = glm::vec3(0);
		X[selected_index].x += Right[0]*valX ;
		float newValue = X[selected_index].y+Up[1]*valY;
		if(newValue>0)
			X[selected_index].y = newValue;
		X[selected_index].z += Right[2]*valX + Up[2]*valY;
	}
	oldX = x;
	oldY = y;

	glutPostRedisplay();
}


void DrawGrid()
{
	glBegin(GL_LINES);
	glColor3f(0.5f, 0.5f, 0.5f);
	for(int i=-GRID_SIZE;i<=GRID_SIZE;i++)
	{
		glVertex3f((float)i,0,(float)-GRID_SIZE);
		glVertex3f((float)i,0,(float)GRID_SIZE);

		glVertex3f((float)-GRID_SIZE,0,(float)i);
		glVertex3f((float)GRID_SIZE,0,(float)i);
	}
	glEnd();
}

void InitGL() {
	startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	currentTime = startTime;

	// get ticks per second
    QueryPerformanceFrequency(&frequency);

    // start timer
    QueryPerformanceCounter(&t1);



	glEnable(GL_DEPTH_TEST);
	int i=0, j=0, count=0;
	int l1=0, l2=0;
	int v = numY+1;
	int u = numX+1;

	indices.resize( numX*numY*2*3);

	X.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);

	V_new.resize(total_points);


	//fill in positions
	for( j=0;j<=numY;j++) {
		for( i=0;i<=numX;i++) {
			X[count++] = glm::vec3( ((float(i)/(u-1)) *2-1)* hsize, size+1, ((float(j)/(v-1) )* size));
		}
	}

	//fill in velocities

	memset(&(V[0].x),0,total_points*sizeof(glm::vec3));

	//fill in indices
	GLushort* id=&indices[0];
	for (i = 0; i < numY; i++) {
		for (j = 0; j < numX; j++) {
			int i0 = i * (numX+1) + j;
			int i1 = i0 + 1;
			int i2 = i0 + (numX+1);
			int i3 = i2 + 1;
			if ((j+i)%2) {
				*id++ = i0; *id++ = i2; *id++ = i1;
				*id++ = i1; *id++ = i2; *id++ = i3;
			} else {
				*id++ = i0; *id++ = i2; *id++ = i3;
				*id++ = i0; *id++ = i3; *id++ = i1;
			}
		}
	}

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glPolygonMode(GL_BACK, GL_LINE);
	glPointSize(5);

	wglSwapIntervalEXT(0);

	//setup springs
	// Horizontal
	for (l1 = 0; l1 < v; l1++)	// v
		for (l2 = 0; l2 < (u - 1); l2++) {
			AddSpring((l1 * u) + l2,(l1 * u) + l2 + 1,KsStruct,KdStruct,STRUCTURAL_SPRING);
		}

	// Vertical
	for (l1 = 0; l1 < (u); l1++)
		for (l2 = 0; l2 < (v - 1); l2++) {
			AddSpring((l2 * u) + l1,((l2 + 1) * u) + l1,KsStruct,KdStruct,STRUCTURAL_SPRING);
		}


	// Shearing Springs
	for (l1 = 0; l1 < (v - 1); l1++)
		for (l2 = 0; l2 < (u - 1); l2++) {
			AddSpring((l1 * u) + l2,((l1 + 1) * u) + l2 + 1,KsShear,KdShear,SHEAR_SPRING);
			AddSpring(((l1 + 1) * u) + l2,(l1 * u) + l2 + 1,KsShear,KdShear,SHEAR_SPRING);
		}


	// Bend Springs
	for (l1 = 0; l1 < (v); l1++) {
		for (l2 = 0; l2 < (u - 2); l2++) {
			AddSpring((l1 * u) + l2,(l1 * u) + l2 + 2,KsBend,KdBend,BEND_SPRING);
		}
		AddSpring((l1 * u) + (u - 3),(l1 * u) + (u - 1),KsBend,KdBend,BEND_SPRING);
	}
	for (l1 = 0; l1 < (u); l1++) {
		for (l2 = 0; l2 < (v - 2); l2++) {
			AddSpring((l2 * u) + l1,((l2 + 2) * u) + l1,KsBend,KdBend,BEND_SPRING);
		}
		AddSpring(((v - 3) * u) + l1,((v - 1) * u) + l1,KsBend,KdBend,BEND_SPRING);
	}

	int total_springs = springs.size();

	M.resize(total_springs);

	M = mass*M;


	C.resize(total_springs );
	C_Dot.resize(total_springs );
	dc_dp.resize(total_springs );
	df_dx.resize(total_springs );
	deltaP2.resize(total_springs );
	memset(&(C[0]),0,total_springs*sizeof(float));
	memset(&(C_Dot[0]),0,total_springs*sizeof(float));
	memset(&(deltaP2[0].x),0,total_springs*sizeof(glm::vec3));

}

void OnReshape(int nw, int nh) {
	glViewport(0,0,nw, nh);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)nw / (GLfloat)nh, 1.f, 100.0f);

	glGetIntegerv(GL_VIEWPORT, viewport);
	glGetDoublev(GL_PROJECTION_MATRIX, P);

	glMatrixMode(GL_MODELVIEW);
}

void OnRender() {
	size_t i=0;
	float newTime = (float) glutGet(GLUT_ELAPSED_TIME);
	frameTime = newTime-currentTime;
	currentTime = newTime;
	//accumulator += frameTime;

	//Using high res. counter
    QueryPerformanceCounter(&t2);
	 // compute and print the elapsed time in millisec
    frameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;
	t1=t2;
	accumulator += frameTimeQP;

	++totalFrames;
	if((newTime-startTime)>1000)
	{
		float elapsedTime = (newTime-startTime);
		fps = (totalFrames/ elapsedTime)*1000 ;
		startTime = newTime;
		totalFrames=0;
	}

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f", fps, frameTime, frameTimeQP);
	glutSetWindowTitle(info);

	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT);
	glLoadIdentity();
	glTranslatef(0,0,dist);
	glRotatef(rX,1,0,0);
	glRotatef(rY,0,1,0);

	glGetDoublev(GL_MODELVIEW_MATRIX, MV);
	viewDir.x = (float)-MV[2];
	viewDir.y = (float)-MV[6];
	viewDir.z = (float)-MV[10];
	Right = glm::cross(viewDir, Up);

	//draw grid
	DrawGrid();

	//draw polygons
	glColor3f(1,1,1);
	glBegin(GL_TRIANGLES);
	for(i=0;i<indices.size();i+=3) {
		glm::vec3 p1 = X[indices[i]];
		glm::vec3 p2 = X[indices[i+1]];
		glm::vec3 p3 = X[indices[i+2]];
		glVertex3f(p1.x,p1.y,p1.z);
		glVertex3f(p2.x,p2.y,p2.z);
		glVertex3f(p3.x,p3.y,p3.z);
	}
	glEnd();

	//draw points

	glBegin(GL_POINTS);
	for(i=0;i<total_points;i++) {
		glm::vec3 p = X[i];
		int is = (i==selected_index);
		glColor3f((float)!is,(float)is,(float)is);
		glVertex3f(p.x,p.y,p.z);
	}
	glEnd();

	glutSwapBuffers();
}

void OnShutdown() {
	X.clear();
	V.clear();
	V_new.clear();
	M.clear();
	df_dx.clear();
	dc_dp.clear();
	F.clear();
	springs.clear();
	indices.clear();
}

void ComputeForces() {
	size_t i=0;

	for(i=0;i<total_points;i++) {
		F[i] = glm::vec3(0);

		//add gravity force
		if(i!=0 && i!=( numX)	)
			F[i] += gravity ;

		//F[i] += DEFAULT_DAMPING*V[i];
	}

	for(i=0;i<springs.size();i++) {
		glm::vec3 p1 = X[springs[i].p1];
		glm::vec3 p2 = X[springs[i].p2];
		glm::vec3 v1 = V[springs[i].p1];
		glm::vec3 v2 = V[springs[i].p2];
		glm::vec3 deltaP = p1-p2;
		glm::vec3 deltaV = v1-v2;
		float dist = glm::length(deltaP);


		//fill in the Jacobian matrix
		//float dist2 = dist*dist;
		//float lo_l  = springs[i].rest_length/dist;
		//K[i] = springs[i].Ks* (-I + (lo_l * (I - glm::dot(deltaP, deltaP)/dist2) ) );

		C[i] = dist-springs[i].rest_length;
		dc_dp[i] = deltaP/dist;
		C_Dot[i] = glm::dot(v1, -dc_dp[i]) + glm::dot(v2, dc_dp[i]);
		deltaP2[i] = glm::vec3(deltaP.x*deltaP.x,deltaP.y*deltaP.y, deltaP.z*deltaP.z);


		float leftTerm = -springs[i].Ks * (dist-springs[i].rest_length);
		float rightTerm = springs[i].Kd * (glm::dot(deltaV, deltaP)/dist);
		glm::vec3 springForce = (leftTerm + rightTerm)*glm::normalize(deltaP);

		if(springs[i].p1 != 0 && springs[i].p1 != numX)
			F[springs[i].p1] += springForce;
		if(springs[i].p2 != 0 && springs[i].p2 != numX )
			F[springs[i].p2] -= springForce;
	}
}
void CalcForceDerivatives() {
	//clear the derivatives
	memset(&(df_dx[0]),0,total_points*sizeof(glm::mat3));
//	memset(&(df_dv[0]),0,total_points*sizeof(glm::mat3));

	size_t i=0;

	glm::mat3 d2C_dp2[2][2]={glm::mat3(1.0f),glm::mat3(1.0f),glm::mat3(1.0f),glm::mat3(1.0f)};

	//#pragma omp parallel for
	for(i=0;i<springs.size();i++) {
		float c1 = C[i];
		d2C_dp2[0][0][0][0] = (-c1*deltaP2[i].x+c1);
		d2C_dp2[0][0][1][1] = (-c1*deltaP2[i].y+c1);
		d2C_dp2[0][0][2][2] = (-c1*deltaP2[i].z+c1);

		d2C_dp2[0][1][0][0] = (c1*deltaP2[i].x-c1);
		d2C_dp2[0][1][1][1] = (c1*deltaP2[i].y-c1);
		d2C_dp2[0][1][2][2] = (c1*deltaP2[i].z-c1);

		d2C_dp2[1][0]  = d2C_dp2[0][1];
		d2C_dp2[1][1]  = d2C_dp2[0][0];

		glm::mat3 dp1  = glm::outerProduct(dc_dp[i], dc_dp[i]);
		glm::mat3 dp2  = glm::outerProduct(dc_dp[i], -dc_dp[i]);
		glm::mat3 dp3  = glm::outerProduct(-dc_dp[i],-dc_dp[i]);

		df_dx[i] += -springs[i].Ks* (dp1 + (d2C_dp2[0][0]*C[i]) ) - springs[i].Kd * (d2C_dp2[0][0]*C_Dot[i]);
		df_dx[i] += -springs[i].Ks* (dp2 + (d2C_dp2[0][1]*C[i]) ) - springs[i].Kd * (d2C_dp2[0][1]*C_Dot[i]);
		df_dx[i] += -springs[i].Ks* (dp2 + (d2C_dp2[1][1]*C[i]) ) - springs[i].Kd * (d2C_dp2[1][1]*C_Dot[i]);

		//df_dv[i] += -springs[i].Kd*dp1;
		//df_dv[i] += -springs[i].Kd*dp2;
		//df_dv[i] += -springs[i].Kd*dp3;
	}

}




void IntegrateImplicit(float deltaTime) {
	float deltaT2 = deltaTime*deltaTime;
	CalcForceDerivatives();

	LargeVector<glm::mat3> A = M - deltaT2* df_dx;
	LargeVector<glm::vec3> b = M*V + deltaTime*F;

	SolveConjugateGradient(A, V_new, b);

	for(size_t i=0;i<total_points;i++) {
		X[i] +=  deltaTime*V_new[i];
		if(X[i].y <0) {
			X[i].y = 0;
		}
		V[i] = V_new[i];
	}
}

void ApplyProvotDynamicInverse() {

	for(size_t i=0;i<springs.size();i++) {
		//check the current lengths of all springs
		glm::vec3 p1 = X[springs[i].p1];
		glm::vec3 p2 = X[springs[i].p2];
		glm::vec3 deltaP = p1-p2;

	 	float dist = glm::length(deltaP);
		if(dist>springs[i].rest_length) {
			dist -= (springs[i].rest_length);
			dist /= 2.0f;
			deltaP = glm::normalize(deltaP);
			deltaP *= dist;
			if(springs[i].p1==0 || springs[i].p1 ==numX) {
				V[springs[i].p2] += deltaP;
			} else if(springs[i].p2==0 || springs[i].p2 ==numX) {
			 	V[springs[i].p1] -= deltaP;
			} else {
				V[springs[i].p1] -= deltaP;
				V[springs[i].p2] += deltaP;
			}
		}
	}
}
void OnIdle() {
	/*
	//Semi-fixed time stepping
	if ( frameTime > 0.0 )
    {
        const float deltaTime = min( frameTime, timeStep );
        StepPhysics(deltaTime );
        frameTime -= deltaTime;
    }
	*/

	//Fixed time stepping + rendering at different fps
	if ( accumulator >= timeStep )
    {
        StepPhysics(timeStep );
        accumulator -= timeStep;
    }
	glutPostRedisplay();
}

void StepPhysics(float dt ) {
	ComputeForces();
		IntegrateImplicit(timeStep);
	ApplyProvotDynamicInverse();
}

void main(int argc, char** argv) {
	atexit(OnShutdown);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("GLUT Cloth Demo [Implicit Euler Integration]");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);

	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);
	glewInit();
	InitGL();

	glutMainLoop();
}
