//This is the simplest possible implementation of the paper
//"Point based animation of elastic, plastic and melting objects" by Matthias Mueller et. al"
//(http://www.matthiasmueller.info/publications/sca04.pdf) and the detailed chapter on
//Meshless Finite Elements, Chapter 7 in Point-Based Graphics,Markus Gross,
//Hanspeter Pfister (eds.), ISBN 0123706041,pp: 341--357, 2007.
//
//This code is under BSD license. If you make some improvements,
//or are using this in your research, do let me know and I would appreciate
//if you acknowledge this in your code.
//
//Controls:
//left click on any empty region to rotate, middle click to zoom
//left click and drag any point to drag it.
//Press '[' or ']' to decrease/increase Poisson Ratio (nu)
//Press ',' or '.' to decrease/increase Young's Modulus (Y)
//Press 'j' to view the Jacobians
//Press 'f' to view the Forces which maybe
//      'x' for displacement,
//      'e' for strains and
//      's' for stresses
//
//Author: Movania Muhammad Mobeen
//        School of Computer Engineering,
//        Nanyang Technological University,
//        Singapore.
//Email : mova0002@e.ntu.edu.sg
//
//Special thanks to Marco Fratarcangeli (http://www.fratarcangeli.net/) for sharing his
//implementation.



#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/glut.h>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>

#include <glm/glm.hpp>

#pragma comment(lib, "glew32.lib")

using namespace std;
const int width = 1024, height = 1024;

int numX = 10, numY=5, numZ=5;
const int total_points = numX*numY*numZ;
int sizeX = 4;
int sizeY = 1;
int sizeZ = 1;
float hsizeX = sizeX/2.0f;
float hsizeY = sizeY/2.0f;
float hsizeZ = sizeZ/2.0f;

bool bShowForces=false, bShowJacobians=false;

float timeStep =  1/60.0f;
float currentTime = 0;
double accumulator = timeStep;
int selected_index = -1;

int what=0;//0-displacements, 1-stresses, 2-strains

struct neighbor { int  j; float w;glm::vec3 rdist, dj; };

glm::mat3		I=glm::mat3(1);			//identity matrix
vector<glm::vec3>			Xi;			//initial positions
vector<glm::vec3>			X;			//current positions
vector<glm::vec3>			U;			//displacements

vector<float>				M;			//masses
vector<glm::mat3>			Minv;		//inverse of moment matrix (Ainv in paper)

vector<glm::vec3>			di;

vector<glm::mat3>			sigma;		//stresses
vector<glm::mat3>			epsilon;	//strains
vector<glm::mat3>			J;			//Jacobian
vector<glm::vec3>			acc0;

vector<vector<neighbor>>	neighbors;	//Neighboring phyxels

vector<float>				r;			//distance to the neighboring phyxel
vector<float>				h;			//support radius of kernel
vector<float>				rho;		//density of phyxel
vector<float>				Vol;		//volume of phyxel
vector<glm::vec3>			V;			//current velocity
vector<glm::vec3>			F;			//total force
vector<bool>				isFixed;

float scFac = 0; //scaling constant


int oldX=0, oldY=0;
float rX=15, rY=0;
int state =1 ;
float dist=-23;
const int GRID_SIZE=10;
float pointSize = 10;
float spacing =  float(sizeX)/(numX+1);							// Spacing of particles


glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);


GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
double frameTimeQP=0;
float frameTime =0 ;

float startTime =0, fps=0;
int totalFrames=0;

int whichIndex = 0;
char info[MAX_PATH]={0};

float nu =	0.4f;				//Poisson ratio
float Y = 300000.0f;			//Young modulus
float density = 10000.f;		//material density
float kv=100, damping=500.0f;	//constant used in Eq. 22 page 5
float d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
float d16 = (1.0f - nu) * d15;
float d17 = nu * d15;
float d18 = Y / 2 / (1.0f + nu);

glm::vec3 D(d16, d17, d18); //Isotropic elasticity matrix D
							//(in the original paper this matrix is called C) in Eq. 4 on page 3


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
		if(winZ==1)
			winZ=0;
		double objX=0, objY=0, objZ=0;
		gluUnProject(window_x,window_y, winZ,  MV,  P, viewport, &objX, &objY, &objZ);
		glm::vec3 pt(objX,objY, objZ);
		int i=0;
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

		V[selected_index].x = 0;
		V[selected_index].y = 0;
		V[selected_index].z = 0;
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


void Print(FILE* fp, glm::mat3 m) {
	for(int j=0;j<3;j++) {
		for(int i=0;i<3;i++) {
			fprintf(fp,"%3.3f ",m[j][i]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");
}

#include <algorithm>
typedef std::pair<int, float> mypair;

struct Cmp {
	bool operator()(const mypair &lhs, const mypair &rhs) {
		return lhs.second < rhs.second;
	}
};

void GetKNearestNeighbors(int index, int k, vector<float>& dis, vector<neighbor>& n) {
	vector<mypair> distances;

	//Get the distances of current point to all other points
	for(int i=0;i<total_points;i++) {
		if(index!=i) {
			mypair m;
			m.first = i;
			m.second= fabs(glm::distance(X[index], X[i]));
			distances.push_back(m);
		}
	}

	//sort the distances
	sort (distances.begin(), distances.end(), Cmp());

	//now take the top k neighbors based on distance
	for(int i=0;i<k;i++) {
		neighbor ne;
		ne.j = distances[i].first;
		dis.push_back(distances[i].second);
		n.push_back(ne);
	}
}
void ComputeRadiusAndSupportRadii(vector<float> &dists, float& r, float& h)
{
	//For this read section 3.2 Initialization on page 4
	//look through all neigbour distances
	//and average the distance. This will give us r
	float avg = 0.f;
	for(size_t i=0;i<dists.size();i++)
		avg += dists[i];
	r = avg / dists.size();

	// compute the support radii H = 3 * R
	h = 3.0f * r;
}
void FillNeighs(int index)
{
	//For this read section 3.2 Initialization on page 4
	//Based on Eq. 9
	float mul_factor = float(315.0f/(64*M_PI*pow(h[index],9.0f)));//35 / (float)(32 * M_PI * pow(h, 7));
	float h2 = h[index]*h[index];
	for(size_t i=0;i<neighbors[index].size();i++)
	{
		neighbor& n = neighbors[index][i];
		n.rdist  = Xi[n.j] - Xi[index];
		float r2 = glm::dot(n.rdist, n.rdist);		//r = sqrt(Xij.x*Xij.x+Xij.y*Xij.y+Xij.z*Xij.z)
													//r^2 = dot(Xij,Xij);
		n.w = mul_factor * pow(h2 - r2, 3 );
	}
}
void ComputeScalingConstant() {
	printf("ComputeScalingFactor: estimating the scaling factor...\n");
	scFac = 0.f;
	int i=0;
	for ( i =0; i<total_points; i++ ) {
		float sum = 0.f;

		vector<neighbor>& pNeighbor = neighbors[i];
		for(size_t j=0;j<pNeighbor.size();j++) {
			sum += pow(r[pNeighbor[j].j], 3) * pNeighbor[j].w;
		}
		scFac += 1.0f / sum;
	}
	// This is the common scaling factor to compute the mass of each phyxel.
	// See last paragraph of Section 3.2 on page 4
	scFac /= total_points;

	printf("Scaling factor: %3.3f\n", scFac);
}
void ComputeMass(float dm, int index)
{
	// See last paragraph of Section 3.2 on page 4
	M[index] = scFac * pow(r[index], 3) * dm;
}
void ComputeDensityAndVolume(int index)
{
	// See last paragraph of Section 3.2 on page 4
	rho[index] = 0.f;
	vector<neighbor>& pNeighbor = neighbors[index];
	for(size_t i=0;i<pNeighbor.size();i++)
		rho[index] += M[pNeighbor[i].j] * pNeighbor[i].w;		// Eq. 10 page 4
	Vol[index] =  M[index] / rho[index];
}

void ComputeInvMomentMatrix(int index)
{
	glm::mat3 A, A_sum, V;
	A =glm::mat3(0);

	vector<neighbor>& pNeighbor = neighbors[index];
	for(size_t i=0;i<pNeighbor.size();i++)
	{
		A_sum = glm::outerProduct(pNeighbor[i].rdist, pNeighbor[i].rdist * pNeighbor[i].w);	// Eq. 14
		A += A_sum;
	}

	if(glm::determinant(A) != 0.0)
		Minv[index] = glm::inverse(A);		// Eq. 14, inverted moment matrix
	else
	{
		// if MomentMatrix is not invertible it means that there are less than 4 neighbors
		// or the neighbor phyxels are coplanar or colinear
		// We should use SVD to extract the inverse but I have left this as is.
		//Read section 3.3 last paragraph
		puts("Warning: Singular matrix!!!");
	}

	di[index]=glm::vec3(0);

	for(size_t i=0;i<pNeighbor.size();i++)
	{
		glm::vec3 Xij_Wij = pNeighbor[i].rdist * pNeighbor[i].w;
		pNeighbor[i].dj = Minv[index] * Xij_Wij;	//Eq. 21 page 5
		di[index] -= pNeighbor[i].dj;				//Eq. 20 page 5
	}
}
void ComputeJacobians()
{
	for(int i=0;i<total_points;i++) {
		vector<neighbor>& pNeighbor = neighbors[i];
		for(size_t j=0;j<pNeighbor.size();j++)
			U[pNeighbor[j].j] = X[pNeighbor[j].j] - Xi[pNeighbor[j].j];

		glm::mat3 B=glm::mat3(0);		// help matrix used to compute the sum in Eq. 15

		// reset du and du_tr
		glm::mat3 du=glm::mat3(0);
		glm::mat3 du_tr=glm::mat3(0);

		for(size_t j=0;j<pNeighbor.size();j++)
		{
			glm::mat3 Bj=glm::mat3(0);
			//Eq. 15 right hand side terms with A_inv
			Bj =glm::outerProduct(U[pNeighbor[j].j] - U[i], pNeighbor[j].rdist * pNeighbor[j].w );
			B += Bj;
		}
		B = glm::transpose(B);

		du = Minv[i] * B;	// Eq. 15 page 4
		du_tr = glm::transpose(du);
		J[i]=glm::mat3(1);
		J[i] += du_tr;		// Eq. 1
	}
}

void ComputeStrainAndStress()
{
	for(int i=0;i<total_points;i++) {
		glm::mat3 Jtr = glm::transpose(J[i]);
		epsilon[i] = (Jtr * J[i]) - I;		// formula 3, Green-Saint Venant non-linear tensor

		glm::mat3& e= epsilon[i];
		glm::mat3& s= sigma[i];

		s[0][0] = D.x*e[0][0]+D.y*e[1][1]+D.y*e[2][2];
		s[1][1] = D.y*e[0][0]+D.x*e[1][1]+D.y*e[2][2];
		s[2][2] = D.y*e[0][0]+D.y*e[1][1]+D.x*e[2][2];

		s[0][1] = D.z*e[0][1];
		s[1][2] = D.z*e[1][2];
		s[2][0] = D.z*e[2][0];

		s[0][2] = s[2][0];
		s[1][0] = s[0][1];
		s[2][1] = s[1][2];
	}
}
void UpdateForces()
{
	//Calculate external force
	for(int i=0;i<total_points;i++) {
		if(!isFixed[i])
			F[i] = glm::vec3(0,gravity.y,0);
		else
			F[i] = glm::vec3(0);

		//Add velocity damping
		F[i] -= V[i] * damping;
	}
	ComputeJacobians();

	ComputeStrainAndStress();


	//Calculate internal force using the stress and Jacobians
	for(int i=0;i<total_points;i++) {
		glm::mat3 F_e, F_v;
		F_e =  -2 * Vol[i] * (J[i] * sigma[i]) ;	// Eq. 18
		glm::vec3 J_u = glm::vec3(J[i][0][0], J[i][0][1],J[i][0][2]);
		glm::vec3 J_v = glm::vec3(J[i][1][0], J[i][1][1],J[i][1][2]);
		glm::vec3 J_w = glm::vec3(J[i][2][0], J[i][2][1],J[i][2][2]);

		glm::vec3 row0 = glm::cross(J_v, J_w);	//Eq. 22
		glm::vec3 row1 = glm::cross(J_w, J_u);	//Eq. 22
		glm::vec3 row2 = glm::cross(J_u, J_v);	//Eq. 22
		glm::mat3 M= glm::transpose(glm::mat3(row0, row1, row2));	//Eq. 22

		F_v = -Vol[i] * kv * (glm::determinant(J[i]) - 1) * M ; //Eq. 22

		vector<neighbor>& pNeighbor = neighbors[i];
		for(size_t j=0;j<pNeighbor.size();j++)
			F[pNeighbor[j].j] += (F_e + F_v) * pNeighbor[j].dj;

		F[i] += (F_e + F_v) * di[i];
	}
}
void InitGL() {
	startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	currentTime = startTime;

	// get ticks per second
	QueryPerformanceFrequency(&frequency);

	// start timer
	QueryPerformanceCounter(&t1);



	glEnable(GL_DEPTH_TEST);
	glEnable(GL_POINT_SMOOTH);
	int i=0, j=0,count=0;

	float ypos = 4.0f;

	Xi.resize(total_points);
	r.resize(total_points);
	h.resize(total_points);
	M.resize(total_points);
	isFixed.resize(total_points);
	J.resize(total_points);
	sigma.resize(total_points);
	epsilon.resize(total_points);
	X.resize(total_points);
	U.resize(total_points);
	V.resize(total_points);
	F.resize(total_points);
	rho.resize(total_points);
	Vol.resize(total_points);
	neighbors.resize(total_points);
	Minv.resize(total_points);
	acc0.resize(total_points);
	di.resize(total_points);

	//fill in X
	for(int k=0;k<numZ;k++) {
		for( j=0;j<numY;j++) {
			for( i=0;i<numX;i++) {
				X[count++] = glm::vec3( ((float(i)/(numX-1)) )*  sizeX,
					((float(j)/(numY-1))*2-1)* hsizeY + ypos,
					((float(k)/(numZ-1))*2-1)* hsizeZ);

			}
		}
	}


	memcpy(&Xi[0].x, &X[0].x, total_points*sizeof(glm::vec3));

	//fill in V
	memset(&(V[0].x),0,total_points*sizeof(glm::vec3));

	for(i=0; i < total_points; ++i) {
		//Fix the phyxels at the X axis edges
		isFixed[i] = (X[i].x==0);//(X[i].x<=-hsizeX || X[i].x>=hsizeX);
	}

	for(i=0; i < total_points; ++i) {
		vector<float> dist;
		GetKNearestNeighbors(i, 10, dist, neighbors[i]);
		ComputeRadiusAndSupportRadii(dist, r[i], h[i]);
	}

	for(i=0;i<total_points;i++) {
		FillNeighs(i);
	}

	ComputeScalingConstant();

	for(i=0; i < total_points; ++i)
		ComputeMass(density, i);

	for(i=0; i < total_points; ++i)
		ComputeDensityAndVolume(i);

	for(i=0; i < total_points; ++i) {
		ComputeInvMomentMatrix(i);
	}
 	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//glPolygonMode(GL_BACK, GL_LINE);
	glPointSize(pointSize );

	wglSwapIntervalEXT(0);
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
	int i=0;
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

	sprintf_s(info, "FPS: %3.2f, Frame time (GLUT): %3.4f ms, (QP): %3.3f ms, Young Mod.: %f, Poisson ratio: %4.4f", fps, frameTime, frameTimeQP, Y, nu);

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

	//draw mesh
	glColor3f(0.75,0.75,0.75);
	glBegin(GL_LINES);
		for(int i=0;i<total_points;i++) {
			for(size_t j=0;j<neighbors[i].size();j++) {
				glVertex3fv(&X[i].x);
				glVertex3fv(&X[neighbors[i][j].j].x);
			}
		}
	glEnd();

	//draw points
	glBegin(GL_POINTS);
	for(i=0;i<total_points;i++) {
		glm::vec3 p = X[i];
		if(i==selected_index)
			glColor3f(0,1,1);
		else
			glColor3f(1,0,0);
		glVertex3f(p.x,p.y,p.z);
	}

	/*
	//Debug code to see the neighbors of a given node
	glColor3f(1,1,0);

	glVertex3fv(&X[whichIndex ].x);
	glColor3f(0,1,1);
	for(int i=0;i<neighbors[whichIndex].size();i++) {
		if(neighbors[whichIndex][i].j !=whichIndex) {
			glVertex3fv(&X[neighbors[whichIndex][i].j].x);
			printf("Index: %d\n",neighbors[whichIndex][i].j);
		}
	}
	 */
	glEnd();


	if(bShowJacobians) {
		glBegin(GL_LINES);
		for(i=0;i<total_points;i++) {
			glm::vec3 p = X[i];
			glm::mat3 j = J[i];

			glColor3f(1,0,0);	glVertex3fv(&p.x);	glVertex3f(p.x+j[0].x,p.y+j[0].y,p.z+j[0].z);
			glColor3f(0,1,0);	glVertex3fv(&p.x);	glVertex3f(p.x+j[1].x,p.y+j[1].y,p.z+j[1].z);
			glColor3f(0,0,1);	glVertex3fv(&p.x);	glVertex3f(p.x+j[2].x,p.y+j[2].y,p.z+j[2].z);
		}
		glEnd();
	}

	if(bShowForces){
		//Visualize displacements
		glColor3f(0,1,0);

		glBegin(GL_LINES);
		for(i=0;i<total_points;i++) {
			glm::vec3 u = glm::vec3(0);
			switch(what) {
				case 0:	u = -U[i]; break;
				case 1: u = glm::normalize(glm::vec3(sigma[i][0][0], sigma[i][0][1], sigma[i][0][2])); break;
				case 2: u = glm::normalize(glm::vec3(epsilon[i][0][0], epsilon[i][0][1], epsilon[i][0][2])) ; break;
			}
			glm::vec3 p = X[i];
			glVertex3fv(&p.x);
			glVertex3f(p.x+u.x,p.y+u.y,p.z+u.z);
		}
		glEnd();
	}

	glutSwapBuffers();
}

void OnShutdown() {

	Xi.clear();
	X.clear();
	U.clear();

	M.clear();
	Minv.clear();

	V.clear();
	F.clear();

	rho.clear();
	Vol.clear();

	di.clear();
	sigma.clear();
	epsilon.clear();
	J.clear();
	acc0.clear();
	neighbors.clear();

	r.clear();
	h.clear();
	isFixed.clear();

}


void IntegrateLeapFrog(float dt) {
	int i=0;
	float dt2 = dt*dt;
	float half_dt2 = dt2*0.5f;

	for(i=0;i<total_points;i++) {

		//X_(i+1) = X_i + V_i*dt + 1/2*a_i*dt^2
		if(!isFixed[i]) {
			acc0[i] = F[i]/M[i];

			X[i] += dt*V[i]+(acc0[i]*half_dt2);

			if(X[i].y <0) {
				X[i].y = 0;
			}
		}
	}
	//Calculate the new acceleration
	UpdateForces();

	//V_(i+1) = V_i + ((a_i+a_(i+1))/2)*dt^2
	for(i=0;i<total_points;i++) {
		if(!isFixed[i])
			V[i] += ((F[i]/M[i] + acc0[i])*half_dt2);
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
	UpdateForces();
	IntegrateLeapFrog(dt);
}

void OnKey(unsigned char k,int , int) {
	switch(k) {
		case 'a':whichIndex--;break;
		case 'd':whichIndex++;break;
		case 'f':bShowForces=!bShowForces;break;
		case 'j':bShowJacobians=!bShowJacobians;break;
		case 'e':what=2;break;
		case 's':what=1;break;
		case 'x':what=0;break;
		case ',':Y-=500;break;
		case '.':Y+=500;break;
		case '[':nu-=0.01f;break;
		case ']':nu+=0.01f;break;
	}
	if(nu>0.49999f)		nu=0.49f;
	if(nu<0)			nu=0;
	if(Y<0.01f)			Y=0.01f;
	if(Y>175000000)		Y=175000000;
	if(whichIndex<0)	whichIndex = 0;

	d15 = Y / (1.0f + nu) / (1.0f - 2 * nu);
	d16 = (1.0f - nu) * d15;
	d17 = nu * d15;
	d18 = Y / 2 / (1.0f + nu);
	D.x=d16;
	D.y=d17;
	D.z=d18;
	whichIndex = (whichIndex%total_points);
	glutPostRedisplay();

}
void main(int argc, char** argv) {
	atexit(OnShutdown);
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("GLUT Meshless FEM Demo using Leap Frog Integration");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);

	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove);

	glutKeyboardFunc(OnKey);
	glewInit();
	InitGL();

	puts("Demo code implementing Meshless FEM");
	puts("By Muhammad Mobeen Movania");
	puts("=====================================");
	puts("Press '[' or ']' to decrease/increase Poisson Ratio (nu)");
	puts("Press ',' or '.' to decrease/increase Young's Modulus (Y)");
	puts("Press 'j' to view the Jacobians");
	puts("Press 'f' to view the Forces\n\tPress 'x' for displacement, 'e' for strains and 's' for stresses");

	glutMainLoop();
}
