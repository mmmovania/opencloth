//A simple soft body implementation using the Co-rotated Linear FEM using Implicit Euler integration 
//based on the outstanding book "Physics based Animation" by Kenny Erleben et. al.
//This code is intended for beginners  so that they may understand what is required 
//to implement a simple soft body simulation using FEM.
//
//This code is under BSD license. If you make some improvements,
//or are using this in your research, do let me know and I would appreciate
//if you acknowledge this in your code.
//
//Controls:
//left click on any empty region to rotate, middle click to zoom 
//left click and drag any point to move it.
//
//Author: Movania Muhammad Mobeen
//        School of Computer Engineering,
//        Nanyang Technological University,
//        Singapore.
//Email : mova0002@e.ntu.edu.sg
//
//Special thanks to:
//Kenny Erleben: From whose book, much of the content in this implementation was understood from.
//Erwin Coumans: From whose implementation of FEM shared here:
//               (http://www.bulletphysics.org/Bullet/phpBB3/viewtopic.php?f=4&t=7241 ), 
//               I was able to understand the theory given in the book Physics based Animation.
//
//Matthias 
//Muller Fischer: Whose contributions (in the form of courses and papers) have been priceless for soft body simulations.

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include <vector>
#include <map>
#include <glm/glm.hpp>


#pragma comment(lib, "glew32.lib")

using namespace std;  
const int width = 1024, height = 1024;

#define EPSILON 0.001f
#define EPS2  EPSILON*EPSILON
size_t total_points=0;
float timeStep =  1/60.0f;
float currentTime = 0;
float accumulator = timeStep;
int selected_index = -1;
 

struct Tetrahedron {
	int indices[4];			//indices
	float volume;			//volume 
	float plastic[6];		//plasticity values
	glm::vec3 e1, e2, e3;	//edges
	glm::mat3 Re;			//Rotational warp of tetrahedron.
	glm::mat3 Ke[4][4];		//Stiffness element matrix
	glm::vec3 B[4];			//Jacobian of shapefunctions; B=SN =[d/dx  0     0 ][wn 0  0]
							//                                  [0    d/dy   0 ][0 wn  0]
							//									[0     0   d/dz][0  0 wn]
							//									[d/dy d/dx   0 ]
							//									[d/dz  0   d/dx]
							//									[0    d/dz d/dy]
};

vector<Tetrahedron> tetrahedra;

float nu = 0.33f;			//Poisson ratio
float Y = 500000.0f;		//Young modulus
float density =1000.0f;
float creep = 0.20f;
float yield = 0.04f;
float mass_damping=1.0f;
float m_max = 0.2f;

float d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
float d16 = (1.0f - nu) * d15;
float d17 = nu * d15;
float d18 = Y / 2 / (1.0f + nu);

glm::vec3 D(d16, d17, d18); //Isotropic elasticity matrix D

int oldX=0, oldY=0;
float rX=15, rY=0;
int state =1 ;
float dist=-2.5f;
const int GRID_SIZE=10;
 

glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);  


GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
float frameTimeQP=0;
float frameTime =0 ;

float startTime =0, fps=0;
int totalFrames=0;


char info[MAX_PATH]={0};

int total_tetrahedra =0;

vector<glm::vec3> Xi;		//Model coordinates
vector<glm::vec3> X;		//Current coordinates
vector<glm::vec3> V;		//Velocity
vector<float> mass;			//Mass matrix
vector<glm::vec3> F;		//Force
vector<bool> IsFixed;		//Fixed point
glm::mat3 I=glm::mat3(1);	//3x3 identity matrix
  
typedef map<int, glm::mat3> matrix_map;
vector<matrix_map> K_row;
vector<matrix_map> A_row;
typedef matrix_map::iterator matrix_iterator;
vector<glm::vec3> F0;
vector<glm::vec3> b;

//For Conjugate Gradient 
vector<glm::vec3> residual;
vector<glm::vec3> prev;
vector<glm::vec3> update;
 
float tiny      = 1e-010f;       // TODO: Should be user controllable
float tolerence = 0.001f;   
int i_max = 20;

bool bUseStiffnessWarping = true;
void UpdateOrientation();
void StepPhysics(float dt);
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
		double objX=0, objY=0, objZ=0;
		gluUnProject(window_x,window_y, winZ,  MV,  P, viewport, &objX, &objY, &objZ);
		glm::vec3 pt(objX,objY, objZ); 
		printf("\nObj [ %3.3f,%3.3f,%3.3f ]",objX, objY, objZ);
		size_t i=0;
		for(i=0;i<total_points;i++) {			 
			if( glm::distance(X[i],pt)<0.01) {
				selected_index = i;
				
				printf("Intersected at %d\n",i);
				printf("Pt [ %3.3f,%3.3f,%3.3f ]\n",X[i].x, X[i].y, X[i].z);
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
		UpdateOrientation();
		glutSetCursor(GLUT_CURSOR_INHERIT);
	}
}  
float GetTetraVolume(glm::vec3 e1, glm::vec3 e2, glm::vec3 e3) {
	return  (glm::dot(e1, glm::cross(e2, e3)) )/ 6.0f;
}
void AddTetrahedron(int i0, int i1, int i2, int i3) {
	Tetrahedron t;

	t.indices[0]=i0;
	t.indices[1]=i1;
	t.indices[2]=i2;
	t.indices[3]=i3;
		
	tetrahedra.push_back(t); 
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
		 
		float delta = 1000/abs(dist);
		float valX = (x - oldX)/delta; 
		float valY = (oldY - y)/delta; 
		if(abs(valX)>abs(valY))
			glutSetCursor(GLUT_CURSOR_LEFT_RIGHT);
		else 
			glutSetCursor(GLUT_CURSOR_UP_DOWN);

		V[selected_index] = glm::vec3(0); 
		X[selected_index].x += Right[0]*valX ;
		float newValue = X[selected_index].y + Up[1]*valY;
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
void CalculateK() {

	for(size_t k=0;k<tetrahedra.size();k++) {
		
		glm::vec3 x0 = Xi[tetrahedra[k].indices[0]];
		glm::vec3 x1 = Xi[tetrahedra[k].indices[1]];
		glm::vec3 x2 = Xi[tetrahedra[k].indices[2]];
		glm::vec3 x3 = Xi[tetrahedra[k].indices[3]];
		
		//For this check page no.: 344-346 of Kenny Erleben's book Physics based Animation
		//Eq. 10.30(a-c)
		glm::vec3 e10 = x1-x0;
		glm::vec3 e20 = x2-x0;
		glm::vec3 e30 = x3-x0;

		tetrahedra[k].e1 = e10;
		tetrahedra[k].e2 = e20;
		tetrahedra[k].e3 = e30;

		tetrahedra[k].volume= GetTetraVolume(e10,e20,e30);
		
		//Eq. 10.32
		glm::mat3 E = glm::mat3(e10.x, e10.y, e10.z,
								e20.x, e20.y, e20.z,
								e30.x, e30.y, e30.z) ;
		float detE = glm::determinant(E);
		float invDetE = 1.0f/detE;	
		
		//Eq. 10.40 (a) & Eq. 10.42 (a)
		//Shape function derivatives wrt x,y,z
		// d/dx N0
		float invE10 = (e20.z*e30.y - e20.y*e30.z)*invDetE;
		float invE20 = (e10.y*e30.z - e10.z*e30.y)*invDetE;
		float invE30 = (e10.z*e20.y - e10.y*e20.z)*invDetE;
		float invE00 = -invE10-invE20-invE30;

		//Eq. 10.40 (b) & Eq. 10.42 (b)
		// d/dy N0
		float invE11 = (e20.x*e30.z - e20.z*e30.x)*invDetE;
		float invE21 = (e10.z*e30.x - e10.x*e30.z)*invDetE;
		float invE31 = (e10.x*e20.z - e10.z*e20.x)*invDetE;
		float invE01 = -invE11-invE21-invE31;

		//Eq. 10.40 (c) & Eq. 10.42 (c)
		// d/dz N0
		float invE12 = (e20.y*e30.x - e20.x*e30.y)*invDetE;
		float invE22 = (e10.x*e30.y - e10.y*e30.x)*invDetE;
		float invE32 = (e10.y*e20.x - e10.x*e20.y)*invDetE;
		float invE02 = -invE12-invE22-invE32;

		glm::vec3 B[4];
		//Eq. 10.43 
		//Bn ~ [bn cn dn]^T
		// bn = d/dx N0 = [ invE00 invE10 invE20 invE30 ]
		// cn = d/dy N0 = [ invE01 invE11 invE21 invE31 ]
		// dn = d/dz N0 = [ invE02 invE12 invE22 invE32 ]
		tetrahedra[k].B[0] = glm::vec3(invE00, invE01, invE02);		
		tetrahedra[k].B[1] = glm::vec3(invE10, invE11, invE12);		
		tetrahedra[k].B[2] = glm::vec3(invE20, invE21, invE22);		
		tetrahedra[k].B[3] = glm::vec3(invE30, invE31, invE32);	 
 
		for(int i=0;i<4;i++) {
			for(int j=0;j<4;j++) {
				glm::mat3& Ke = tetrahedra[k].Ke[i][j];
				float d19 = tetrahedra[k].B[i].x;
				float d20 = tetrahedra[k].B[i].y;
				float d21 = tetrahedra[k].B[i].z;
				float d22 = tetrahedra[k].B[j].x;
				float d23 = tetrahedra[k].B[j].y;
				float d24 = tetrahedra[k].B[j].z;
				Ke[0][0]= d16 * d19 * d22 + d18 * (d20 * d23 + d21 * d24);
				Ke[0][1]= d17 * d19 * d23 + d18 * (d20 * d22);
				Ke[0][2]= d17 * d19 * d24 + d18 * (d21 * d22);

				Ke[1][0]= d17 * d20 * d22 + d18 * (d19 * d23);
				Ke[1][1]= d16 * d20 * d23 + d18 * (d19 * d22 + d21 * d24);
				Ke[1][2]= d17 * d20 * d24 + d18 * (d21 * d23);

				Ke[2][0]= d17 * d21 * d22 + d18 * (d19 * d24);
				Ke[2][1]= d17 * d21 * d23 + d18 * (d20 * d24);
				Ke[2][2]= d16 * d21 * d24 + d18 * (d20 * d23 + d19 * d22);

				Ke *= tetrahedra[k].volume;				 
			}
		}
 	}
}

 
void GenerateBlocks(size_t xdim, size_t ydim, size_t zdim, float width, float height, float depth) {
	total_points = (xdim+1)*(ydim+1)*(zdim+1);
	X.resize(total_points);
	Xi.resize(total_points);
	IsFixed.resize(total_points);
	int ind=0;
	float hzdim = zdim/2.0f;
	for(size_t x = 0; x <= xdim; ++x) {
		for (unsigned int y = 0; y <= ydim; ++y) {
			for (unsigned int z = 0; z <= zdim; ++z) {			  
				X[ind] = glm::vec3(width*x,height*z,depth*y);            
				Xi[ind] = X[ind];
				 
				//Make the first few points fixed
				if(Xi[ind].x<0.01)
					IsFixed[ind]=true;

				ind++;
			}
		}
	}
	//offset the tetrahedral mesh by 0.5 units on y axis
	//and 0.5 of the depth in z axis
	for(size_t i=0;i<total_points;i++) {
		X[i].y += 0.5;		
		X[i].z -= hzdim*depth; 
	} 

	for (size_t i = 0; i < xdim; ++i) {
		for (size_t j = 0; j < ydim; ++j) {
			for (size_t k = 0; k < zdim; ++k) {
				int p0 = (i * (ydim + 1) + j) * (zdim + 1) + k;
				int p1 = p0 + 1;
				int p3 = ((i + 1) * (ydim + 1) + j) * (zdim + 1) + k;
				int p2 = p3 + 1;
				int p7 = ((i + 1) * (ydim + 1) + (j + 1)) * (zdim + 1) + k;
				int p6 = p7 + 1;
				int p4 = (i * (ydim + 1) + (j + 1)) * (zdim + 1) + k;
				int p5 = p4 + 1;
				// Ensure that neighboring tetras are sharing faces
				if ((i + j + k) % 2 == 1) {
					AddTetrahedron(p1,p2,p6,p3);
					AddTetrahedron(p3,p6,p4,p7);
					AddTetrahedron(p1,p4,p6,p5);
					AddTetrahedron(p1,p3,p4,p0);
					AddTetrahedron(p1,p6,p4,p3); 
				} else {
					AddTetrahedron(p2,p0,p5,p1);
					AddTetrahedron(p2,p7,p0,p3);
					AddTetrahedron(p2,p5,p7,p6);
					AddTetrahedron(p0,p7,p5,p4);
					AddTetrahedron(p2,p0,p7,p5); 
				}
				total_tetrahedra+=5;
			}
		}
	}

}
void RecalcMassMatrix() {
	//This is a lumped mass matrix
	//Based on Eq. 10.106 and pseudocode in Fig. 10.9 on page 358
	for(size_t i=0;i<total_points;i++) {
		if(IsFixed[i])
			mass[i] = std::numeric_limits<float>::max();
		else
			mass[i] = 1.0f/total_points;
	}

	for(int i=0;i<total_tetrahedra;i++) {
		float m = (density*tetrahedra[i].volume)* 0.25f;				 
		mass[tetrahedra[i].indices[0]] += m;
		mass[tetrahedra[i].indices[1]] += m;
		mass[tetrahedra[i].indices[2]] += m;
		mass[tetrahedra[i].indices[3]] += m;
	}	 
}
 
void InitializePlastic() {
	for(size_t i=0;i<tetrahedra.size();i++) {
		for(int j=0;j<6;j++) 
			tetrahedra[i].plastic[j]=0;		
	} 
}
 
void ClearStiffnessAssembly() {	 
	for(size_t k=0;k<total_points;k++) {
		F0[k].x=0.0f;
		F0[k].y=0.0f;
		F0[k].z=0.0f;
		
		for (matrix_iterator Kij = K_row[k].begin() ; Kij != K_row[k].end(); ++Kij )
			Kij->second=glm::mat3(0);
	}
}
void InitGL() {  
	 
	GenerateBlocks(10,3,3,0.1f,0.1f,0.1f);
	total_tetrahedra = tetrahedra.size();	
 
	total_points = X.size();
	mass.resize(total_points);
 
	//copy positions to buffer 
	A_row.resize(total_points);
	K_row.resize(total_points);
	b.resize(total_points);
	V.resize(total_points);
	F.resize(total_points); 
	F0.resize(total_points); 
	residual.resize(total_points);
	update.resize(total_points);
	prev.resize(total_points);

	//fill in V
	memset(&(V[0].x),0,total_points*sizeof(glm::vec3));

	startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	currentTime = startTime;

	// get ticks per second
	QueryPerformanceFrequency(&frequency);

	// start timer
	QueryPerformanceCounter(&t1);
	glEnable(GL_DEPTH_TEST);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glPointSize(5);
	wglSwapIntervalEXT(0);

	CalculateK();
	ClearStiffnessAssembly();
	RecalcMassMatrix();
	InitializePlastic();
}

void OnReshape(int nw, int nh) {
	glViewport(0,0,nw, nh);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)nw / (GLfloat)nh, 0.1f, 100.0f);

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
	frameTimeQP = (t2.QuadPart - t1.QuadPart) * 1000.0f / frequency.QuadPart;
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

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f, Stiffness Warp: %s", fps, frameTime, frameTimeQP, bUseStiffnessWarping?"On":"Off");
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
 
	glColor3f(0.75,0.75,0.75);	
	glBegin(GL_LINES);

	for(int i=0;i<total_tetrahedra;i++) {		
		int i0 = tetrahedra[i].indices[0];
		int i1 = tetrahedra[i].indices[1];
		int i2 = tetrahedra[i].indices[2];
		int i3 = tetrahedra[i].indices[3];
		glm::vec3 p1 = X[i0];
		glm::vec3 p2 = X[i1];		
		glm::vec3 p3 = X[i2];		
		glm::vec3 p4 = X[i3];

		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p1.x, p1.y, p1.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p4.x, p4.y, p4.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p2.x, p2.y, p2.z);
		glVertex3f(p1.x, p1.y, p1.z);		glVertex3f(p3.x, p3.y, p3.z);

		glVertex3f(p2.x, p2.y, p2.z);		glVertex3f(p3.x, p3.y, p3.z);
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
	Xi.clear();
	V.clear();
	mass.clear();
	F.clear();
	IsFixed.clear();
	tetrahedra.clear();
	K_row.clear();
	A_row.clear();
	F0.clear();
	b.clear();
	residual.clear();
	prev.clear();
	update.clear();
}

void ComputeForces() {
	size_t i=0; 
	for(i=0;i<total_points;i++) {
		F[i] = glm::vec3(0);

		//add gravity force only for non-fixed points
		F[i] += gravity*mass[i];
	}
}


glm::mat3 ortho_normalize(glm::mat3 A) {
	glm::vec3 row0(A[0][0],A[0][1],A[0][2]);
	glm::vec3 row1(A[1][0],A[1][1],A[1][2]);
	glm::vec3 row2(A[2][0],A[2][1],A[2][2]);

	float L0 = glm::length(row0);
	if(L0) 
		row0 /= L0;

	row1 -=  row0 * glm::dot(row0 , row1);
	float L1 = glm::length(row1);
	if(L1) 
		row1 /= L1;

	row2 = glm::cross( row0 , row1);

	return glm::mat3(row0,
					 row1,
					 row2);
}

void UpdateOrientation() {
	for(int k=0;k<total_tetrahedra;k++) {
		//Based on description on page 362-364 
		float div6V = 1.0f / tetrahedra[k].volume*6.0f;

		int i0 = tetrahedra[k].indices[0];
		int i1 = tetrahedra[k].indices[1];
		int i2 = tetrahedra[k].indices[2];
		int i3 = tetrahedra[k].indices[3];

		glm::vec3 p0 = X[i0];
		glm::vec3 p1 = X[i1];
		glm::vec3 p2 = X[i2];
		glm::vec3 p3 = X[i3];

		glm::vec3 e1 = p1-p0;
		glm::vec3 e2 = p2-p0;
		glm::vec3 e3 = p3-p0;

		//Eq. 10.129 on page 363 & Eq. 10.131 page 364
		//n1,n2,n3 approximate Einv
		glm::vec3 n1 = glm::cross(e2,e3) * div6V;
		glm::vec3 n2 = glm::cross(e3,e1) * div6V;
		glm::vec3 n3 = glm::cross(e1,e2) * div6V;
		
		//Now get the rotation matrix from the initial undeformed (model/material coordinates)
		//We get the precomputed edge values
		e1 = tetrahedra[k].e1;
		e2 = tetrahedra[k].e2;
		e3 = tetrahedra[k].e3;

		//Based on Eq. 10.133		
		
		tetrahedra[k].Re[0][0] = e1.x * n1.x + e2.x * n2.x + e3.x * n3.x;  
		tetrahedra[k].Re[0][1] = e1.x * n1.y + e2.x * n2.y + e3.x * n3.y;   
		tetrahedra[k].Re[0][2] = e1.x * n1.z + e2.x * n2.z + e3.x * n3.z;

        tetrahedra[k].Re[1][0] = e1.y * n1.x + e2.y * n2.x + e3.y * n3.x;  
		tetrahedra[k].Re[1][1] = e1.y * n1.y + e2.y * n2.y + e3.y * n3.y;   
		tetrahedra[k].Re[1][2] = e1.y * n1.z + e2.y * n2.z + e3.y * n3.z;

        tetrahedra[k].Re[2][0] = e1.z * n1.x + e2.z * n2.x + e3.z * n3.x;  
		tetrahedra[k].Re[2][1] = e1.z * n1.y + e2.z * n2.y + e3.z * n3.y;  
		tetrahedra[k].Re[2][2] = e1.z * n1.z + e2.z * n2.z + e3.z * n3.z;
		
		tetrahedra[k].Re = ortho_normalize(tetrahedra[k].Re);
		
	}
}
void ResetOrientation() {	
	for(int k=0;k<total_tetrahedra;k++) {
		tetrahedra[k].Re = I;
	}
}
 
void StiffnessAssembly() { 
	
	for(int k=0;k<total_tetrahedra;k++) {
		glm::mat3 Re = tetrahedra[k].Re;
		glm::mat3 ReT = glm::transpose(Re);
 

		for (int i = 0; i < 4; ++i) {
			//Based on pseudocode given in Fig. 10.11 on page 361
			glm::vec3 f=glm::vec3(0.0f,0.0f,0.0f);
			for (int j = 0; j < 4; ++j) {
				glm::mat3 tmpKe = tetrahedra[k].Ke[i][j];
				glm::vec3 x0 = Xi[tetrahedra[k].indices[j]];
				glm::vec3 prod = glm::vec3(tmpKe[0][0]*x0.x+ tmpKe[0][1]*x0.y+tmpKe[0][2]*x0.z, //tmpKe*x0;
										   tmpKe[1][0]*x0.x+ tmpKe[1][1]*x0.y+tmpKe[1][2]*x0.z,
										   tmpKe[2][0]*x0.x+ tmpKe[2][1]*x0.y+tmpKe[2][2]*x0.z);
				f += prod;				   
				if (j >= i) {
					//Based on pseudocode given in Fig. 10.12 on page 361
					glm::mat3 tmp = Re*tmpKe*ReT;
					int index = tetrahedra[k].indices[i]; 		
					 
					K_row[index][tetrahedra[k].indices[j]]+=(tmp);
					 
					if (j > i) {
						index = tetrahedra[k].indices[j];
						K_row[index][tetrahedra[k].indices[i]]+=(glm::transpose(tmp));
					}
				}

			}
			int idx = tetrahedra[k].indices[i];
			F0[idx] -= Re*f;		
		}  	
	} 
}

void AddPlasticityForce(float dt) {
	for(int k=0;k<total_tetrahedra;k++) {
		float e_total[6];
		float e_elastic[6];
		for(int i=0;i<6;++i)
			e_elastic[i] = e_total[i] = 0;

		//--- Compute total strain: e_total  = Be (Re^{-1} x - x0)
		for(unsigned int j=0;j<4;++j) {

			glm::vec3 x_j  =  X[tetrahedra[k].indices[j]];
			glm::vec3 x0_j = Xi[tetrahedra[k].indices[j]];
			glm::mat3 ReT  = glm::transpose(tetrahedra[k].Re);
			glm::vec3 prod = glm::vec3(ReT[0][0]*x_j.x+ ReT[0][1]*x_j.y+ReT[0][2]*x_j.z, //tmpKe*x0;
									   ReT[1][0]*x_j.x+ ReT[1][1]*x_j.y+ReT[1][2]*x_j.z,
									   ReT[2][0]*x_j.x+ ReT[2][1]*x_j.y+ReT[2][2]*x_j.z);
				
			glm::vec3 tmp = prod - x0_j;

			//B contains Jacobian of shape funcs. B=SN
			float bj = tetrahedra[k].B[j].x;
			float cj = tetrahedra[k].B[j].y;
			float dj = tetrahedra[k].B[j].z;

			e_total[0] += bj*tmp.x;
			e_total[1] +=            cj*tmp.y;
			e_total[2] +=                       dj*tmp.z;
			e_total[3] += cj*tmp.x + bj*tmp.y;
			e_total[4] += dj*tmp.x            + bj*tmp.z;
			e_total[5] +=            dj*tmp.y + cj*tmp.z;
		}

		//--- Compute elastic strain
		for(int i=0;i<6;++i)
			e_elastic[i] = e_total[i] - tetrahedra[k].plastic[i];

		//--- if elastic strain exceeds c_yield then it is added to plastic strain by c_creep
		float norm_elastic = 0;
		for(int i=0;i<6;++i)
			norm_elastic += e_elastic[i]*e_elastic[i];
		norm_elastic = sqrt(norm_elastic);
		if(norm_elastic > yield) {
			float amount = dt*min(creep,(1.0f/dt));  //--- make sure creep do not exceed 1/dt
			for(int i=0;i<6;++i)
				tetrahedra[k].plastic[i] += amount*e_elastic[i];
		}

		//--- if plastic strain exceeds c_max then it is clamped to maximum magnitude
		float norm_plastic = 0;
		for(int i=0;i<6;++i)
			norm_plastic += tetrahedra[k].plastic[i]* tetrahedra[k].plastic[i];
		norm_plastic = sqrt(norm_plastic);

		if(norm_plastic > m_max) { 
			float scale = m_max/norm_plastic;
			for(int i=0;i<6;++i)
				tetrahedra[k].plastic[i] *= scale;
		}

		for(size_t n=0;n<4;++n) {
			float* e_plastic = tetrahedra[k].plastic;
			//bn, cn and dn are the shape function derivative wrt. x,y and z axis
			//These were calculated in CalculateK function

			//Eq. 10.140(a) & (b) on page 365
			float bn = tetrahedra[k].B[n].x;
			float cn = tetrahedra[k].B[n].y;
			float dn = tetrahedra[k].B[n].z;
			float D0 = D.x;
			float D1 = D.y;
			float D2 = D.z;
			glm::vec3 f  = glm::vec3(0);

			float  bnD0 = bn*D0;
			float  bnD1 = bn*D1;
			float  bnD2 = bn*D2;
			float  cnD0 = cn*D0;
			float  cnD1 = cn*D1;
			float  cnD2 = cn*D2;
			float  dnD0 = dn*D0;
			float  dnD1 = dn*D1;
			float  dnD2 = dn*D2;
			
			//Eq. 10.141 on page 365
			f.x = bnD0*e_plastic[0] + bnD1*e_plastic[1] + bnD1*e_plastic[2] + cnD2*e_plastic[3] + dnD2*e_plastic[4];
			f.y = cnD1*e_plastic[0] + cnD0*e_plastic[1] + cnD1*e_plastic[2] + bnD2*e_plastic[3] +                  + dnD2*e_plastic[5];
			f.z = dnD1*e_plastic[0] + dnD1*e_plastic[1] + dnD0*e_plastic[2] +                    bnD2*e_plastic[4] + cnD2*e_plastic[5];
			
			f *= tetrahedra[k].volume;
			int idx = tetrahedra[k].indices[n];
			F[idx] += tetrahedra[k].Re*f;
		}
	}
}
void DynamicsAssembly(float dt) {
	float dt2 = dt*dt;
	
	for(size_t k=0;k<total_points;k++) {

		float m_i = mass[k];
		b[k]=glm::vec3(0.0,0.0,0.0);
		//printf("Idx: %3d\n",k);
		
		matrix_map tmp = K_row[k];
		matrix_iterator Kbegin = tmp.begin();
        matrix_iterator Kend   = tmp.end();
		for (matrix_iterator K = Kbegin; K != Kend;++K)
		{
            unsigned int j  = K->first;
			glm::mat3 K_ij  = K->second;
			glm::vec3 x_j   = X[j];	
			glm::mat3& A_ij = A_row[k][j];
 
			A_ij = K_ij * (dt2); 
			glm::vec3 prod = glm::vec3(	K_ij[0][0]*x_j.x + K_ij[0][1]*x_j.y + K_ij[0][2]*x_j.z, 
										K_ij[1][0]*x_j.x + K_ij[1][1]*x_j.y + K_ij[1][2]*x_j.z,
										K_ij[2][0]*x_j.x + K_ij[2][1]*x_j.y + K_ij[2][2]*x_j.z);

            b[k] -= prod;//K_ij * x_j;
			 
            if (k == j)
            {
              float c_i = mass_damping*m_i;
              float tmp = m_i + dt*c_i;
              A_ij[0][0] += tmp; 
			  A_ij[1][1] += tmp;  
			  A_ij[2][2] += tmp;
			}
		}
	  
		b[k] -= F0[k];
		b[k] += F[k];
		b[k] *= dt;
		b[k] += V[k]*m_i;
	} 
}

 



void ConjugateGradientSolver(float dt) {

	
	for(size_t k=0;k<total_points;k++) {
		if(IsFixed[k])
			continue;
		residual[k] = b[k];
 
		matrix_iterator Abegin = A_row[k].begin();
        matrix_iterator Aend   = A_row[k].end();
		for (matrix_iterator A = Abegin; A != Aend;++A)
		{
            unsigned int j   = A->first;
			glm::mat3& A_ij  = A->second;
			float v_jx = V[j].x;	
			float v_jy = V[j].y;
			float v_jz = V[j].z;
			glm::vec3 prod = glm::vec3(	A_ij[0][0] * v_jx+A_ij[0][1] * v_jy+A_ij[0][2] * v_jz, //A_ij * prev[j]
										A_ij[1][0] * v_jx+A_ij[1][1] * v_jy+A_ij[1][2] * v_jz,			
										A_ij[2][0] * v_jx+A_ij[2][1] * v_jy+A_ij[2][2] * v_jz);
			residual[k] -= prod;//  A_ij * v_j;
			
		}
		prev[k]=residual[k];
	}
	
	for(int i=0;i<i_max;i++) {
		float d =0;
		float d2=0;
		
	 	for(size_t k=0;k<total_points;k++) {

			if(IsFixed[k])
				continue;

			update[k].x=0;update[k].y=0;update[k].z=0;
			 
			matrix_iterator Abegin = A_row[k].begin();
			matrix_iterator Aend   = A_row[k].end();
			for (matrix_iterator A = Abegin; A != Aend;++A) {
				unsigned int j   = A->first;
				glm::mat3& A_ij  = A->second;
				float prevx = prev[j].x;
				float prevy = prev[j].y;
				float prevz = prev[j].z;
				glm::vec3 prod = glm::vec3(	A_ij[0][0] * prevx+A_ij[0][1] * prevy+A_ij[0][2] * prevz, //A_ij * prev[j]
											A_ij[1][0] * prevx+A_ij[1][1] * prevy+A_ij[1][2] * prevz,			
											A_ij[2][0] * prevx+A_ij[2][1] * prevy+A_ij[2][2] * prevz);
				update[k] += prod;//A_ij*prev[j];
				 
			}
			d += glm::dot(residual[k],residual[k]);
			d2 += glm::dot(prev[k],update[k]);
		} 
		
		if(fabs(d2)<tiny)
			d2 = tiny;

		float d3 = d/d2;
		float d1 = 0;

		
		for(size_t k=0;k<total_points;k++) {
			if(IsFixed[k])
				continue;

			V[k] += prev[k]* d3;
			residual[k] -= update[k]*d3;
			d1 += glm::dot(residual[k],residual[k]);
		}
		
		if(i >= i_max && d1 < tolerence)
			break;

		if(fabs(d)<tiny)
			d = tiny;

		float d4 = d1/d;
		
		for(size_t k=0;k<total_points;k++) {
			if(IsFixed[k])
				continue;
			prev[k] = residual[k] + prev[k]*d4;
		}		
	}	
}

void UpdatePosition(float dt) {
	for(size_t k=0;k<total_points;k++) {
		if(IsFixed[k])
			continue;
		X[k] += float(dt)*V[k];
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
void GroundCollision()  
{
	for(size_t i=0;i<total_points;i++) {	
		if(X[i].y<0) //collision with ground
			X[i].y=0;
	}
}

void StepPhysics(float dt ) {

	ComputeForces();

	ClearStiffnessAssembly();	

	if(bUseStiffnessWarping)
		UpdateOrientation();
	else
		ResetOrientation();

	StiffnessAssembly();
 
	AddPlasticityForce(dt);
 
	DynamicsAssembly(dt);
 
	ConjugateGradientSolver(dt);
 
	UpdatePosition(dt);

	GroundCollision(); 
}

void OnKey(unsigned char key, int, int) {
	switch(key) {
		case ' ': bUseStiffnessWarping=!bUseStiffnessWarping;	break;
	}
	printf("Stiffness Warping %s\n",bUseStiffnessWarping?"On":"Off");
	glutPostRedisplay();
}

void main(int argc, char** argv) {
	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("GLUT OpenCloth - Co-rotated Linear FEM with Stiffness Warping");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);

	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);	

	glutKeyboardFunc(OnKey);
	glutCloseFunc(OnShutdown);

	glewInit();
	InitGL();

	puts("Press ' ' to toggle stiffness warping on/off\n");

	glutMainLoop();		
}
