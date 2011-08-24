
//A simple cloth using Verlet integration based on the SIGGRAPH course notes
//"Realtime Physics" http://www.matthiasmueller.info/realtimephysics/coursenotes.pdf
//using GLUT,GLEW and GLM libraries. This code is intended for beginners 
//so that they may understand what is required to do Verlet integration
//based cloth simulation. There are two modes in this code:
// 1) CPU mode whereby the Verlet integration is carried out on the CPU as in an 
//    earlier code.
// 2) GPU mode whereby the Verlet integration is carried out on the GPU using the GPGPU
//    approach with GLSL as given in this (http://wwwcg.in.tum.de/Research/data/Publications/simpra05.pdf) 
//    paper. We adopt the point centric approach. Refer to the original paper for details.
//
//This code is under BSD license. If you make some improvements,
//or are using this in your research, do let me know and I would appreciate
//if you acknowledge this in your code.
//
//Controls:
//   Spacebar to toggle btw the two modes i.e. CPU and GPU
//   left click on any empty region to rotate, middle click to zoom 
//   left click and drag any point to drag it.
//
//Author: Movania Muhammad Mobeen
//        School of Computer Engineering,
//        Nanyang Technological University,
//        Singapore.
//Email : mova0002@e.ntu.edu.sg
//Last Modified: 24th August 2011.

#include <GL/glew.h>
#include <GL/wglew.h>
#include <GL/freeglut.h>
#include "GLSLShader.h"
#include <vector>
#include <glm/glm.hpp>
#include <cassert>
using namespace std;  
#pragma comment(lib, "glew32.lib")

#define CHECK_GL_ERRORS assert(glGetError()==GL_NO_ERROR);

const int width = 1024, height = 1024;

int numX = 31, numY=31;
const size_t total_points = (numX+1)*(numY+1);
int sizeX = 4, 
	sizeY = 4;
float hsize = sizeX/2.0f;

const int NUM_ITER = 5;
int selected_index = -1;

enum Mode {CPU, GPU};

struct Spring {
	int p1, p2;
	float rest_length;
	float Ks, Kd;
	int type;
};

Mode current_mode = GPU;
vector<GLushort> indices;
vector<Spring> springs;

vector<glm::vec4> X;
vector<glm::vec4> X_last;
vector<glm::vec3> F;

int oldX=0, oldY=0;
float rX=15, rY=0;
int state =1 ;
float dist=-23;
const int GRID_SIZE=10;

const int STRUCTURAL_SPRING = 0;
const int SHEAR_SPRING = 1;
const int BEND_SPRING = 2;
int spring_count=0;

char info[MAX_PATH]={0};

const float DEFAULT_DAMPING =  -0.125f;
float	KsStruct = 50.75f,KdStruct = -0.25f; 
float	KsShear = 50.75f,KdShear = -0.25f;
float	KsBend = 50.95f,KdBend = -0.25f;
glm::vec3 gravity=glm::vec3(0.0f,-9.81f,0.0f);  
float mass = 10.0f/total_points;

float timeStep =  1.0f/60.0f;
float currentTime = 0;
double accumulator = timeStep;

GLint viewport[4];
GLdouble MV[16];
GLdouble P[16];

glm::vec3 Up=glm::vec3(0,1,0), Right, viewDir;

LARGE_INTEGER frequency;        // ticks per second
LARGE_INTEGER t1, t2;           // ticks
double frameTimeQP=0;
float frameTime =0 ;
int texture_size_x=0;
int texture_size_y=0;

float startTime =0, fps=0 ;
int totalFrames=0;


GLSLShader verletShader, renderShader;
GLuint fboID[2];
GLuint attachID[4];
int readID=0, writeID = 1;
GLuint vboID;
GLenum mrt[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};


float* data[2];
size_t i=0;
GLfloat vRed[4]={1,0,0,1};
GLfloat vWhite[4]={1,1,1,1};

glm::vec3 vec3(glm::vec4 v) {
	return glm::vec3(v.x, v.y, v.z);
}

void InitFBO() { 
	data[0] = &X[0].x;
	data[1] = &X_last[0].x;
	glGenTextures(4, attachID);
	glGenFramebuffers(2, fboID);
	
	for(int j=0;j<2;j++) {
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fboID[j]);	
		printf("FBO %d has textures: ", j);
		for(int i=0;i<2;i++) {
			glBindTexture(GL_TEXTURE_2D, attachID[i+2*j]);
			
			glPixelStorei(GL_UNPACK_ALIGNMENT,1);	//not needed here since our data is 32 bits aligned but added here 
													//for consistency
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
			glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, texture_size_x, texture_size_y, 0, GL_RGBA, GL_FLOAT, data[i]); // NULL = Empty texture
			
			glFramebufferTexture2D(GL_DRAW_FRAMEBUFFER, mrt[i],	GL_TEXTURE_2D, attachID[i+2*j], 0);			 
			printf(" %d ", i+ 2*j);
		}
		printf("attached\n");
	}
	 
	 
	GLenum status = glCheckFramebufferStatus(GL_DRAW_FRAMEBUFFER);
	if (status == GL_FRAMEBUFFER_COMPLETE ) {
		printf("FBO setup succeeded.");
	} else {
		printf("Problem with FBO setup.");
	}
	
	glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
	CHECK_GL_ERRORS
}

void StepPhysics(float dt );
void AddSpring(int a, int b, float ks, float kd, int type) {
	Spring spring;
	spring.p1=a;
	spring.p2=b;
	spring.Ks=ks;
	spring.Kd=kd;
	spring.type = type;
	glm::vec3 deltaP = vec3(X[a]-X[b]);
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
		if(current_mode == CPU) {
			for(i=0;i<total_points;i++) {			 
				if( glm::distance(vec3(X[i]),pt)<0.1) {
					selected_index = i;
					printf("Intersected at %d\n",i);
					break;
				}
			}
		} else {
			glBindBuffer(GL_ARRAY_BUFFER, vboID);
			glm::vec4* ptr = (glm::vec4*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
			// if the pointer is valid(mapped), update VBO
			if(ptr) {
				for(i=0;i<total_points;i++) {
					if( glm::distance(glm::vec3(ptr[i].x, ptr[i].y, ptr[i].z),pt)<0.1) {
						selected_index = i;
						printf("Intersected at %d\n",i);
						break;
					}
				}
				glUnmapBuffer(GL_ARRAY_BUFFER); // unmap it after use
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

		if(current_mode==CPU) {
			X[selected_index].x += Right[0]*valX ;
			float newValue = X[selected_index].y+Up[1]*valY;
			if(newValue>0)
				X[selected_index].y = newValue;
			X[selected_index].z += Right[2]*valX + Up[2]*valY;		
			X_last[selected_index] = X[selected_index];
		} else {
			
			glm::vec4* ptr = (glm::vec4*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
			glm::vec4 oldVal = ptr[selected_index];
			glUnmapBuffer(GL_ARRAY_BUFFER); // unmap it after use

			glm::vec4 newVal;
			newVal.w =1;
			// if the pointer is valid(mapped), update VBO
			if(ptr) {
				// modify buffer data				
				oldVal.x += Right[0]*valX ;
				
				float newValue = oldVal.y+Up[1]*valY;
				if(newValue>0)
					oldVal.y = newValue;
				oldVal.z += Right[2]*valX + Up[2]*valY;					
				newVal=oldVal;				
			}	
			
			int xoff = selected_index%texture_size_x;
			int yoff = selected_index/texture_size_x;
			
			glBindTexture(GL_TEXTURE_2D, attachID[2*readID]);
			glTexSubImage2D(GL_TEXTURE_2D,0,xoff,yoff,1,1, GL_RGBA, GL_FLOAT, &newVal[0]);

			glBindTexture(GL_TEXTURE_2D, attachID[2*readID+1]);
			glTexSubImage2D(GL_TEXTURE_2D,0,xoff,yoff,1,1, GL_RGBA, GL_FLOAT, &newVal[0]);
			
		}
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

void InitVBO(){
	const int size = texture_size_x*texture_size_y*4*sizeof(float);
	glGenBuffers(1, &vboID);
	glBindBuffer(GL_ARRAY_BUFFER, vboID);
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW); 
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}


void InitGL() { 
	texture_size_x =  numX+1;
	texture_size_y =  numY+1;
	CHECK_GL_ERRORS

	renderShader.LoadFromFile(GL_VERTEX_SHADER, "shaders/render.vs");
	renderShader.LoadFromFile(GL_FRAGMENT_SHADER, "shaders/render.fs");
	renderShader.CreateAndLinkProgram();
	renderShader.Use();
		renderShader.AddUniform("color"); 
		renderShader.AddUniform("selected_index"); 
	renderShader.UnUse(); 
	CHECK_GL_ERRORS

	verletShader.LoadFromFile(GL_VERTEX_SHADER, "shaders/verlet.vs");
	verletShader.LoadFromFile(GL_FRAGMENT_SHADER, "shaders/verlet.fs");
	verletShader.CreateAndLinkProgram();
	verletShader.Use();		 
		verletShader.AddUniform("X");				glUniform1i(verletShader("X"),0);		//current position sampler
		verletShader.AddUniform("X_last");			glUniform1i(verletShader("X_last"),1);	//previous position sampler
		verletShader.AddUniform("DEFAULT_DAMPING");	glUniform1f(verletShader("DEFAULT_DAMPING"),DEFAULT_DAMPING);
		verletShader.AddUniform("inv_mass");		glUniform1f(verletShader("inv_mass"),mass);
		verletShader.AddUniform("dt");				glUniform1f(verletShader("dt"),timeStep);
		verletShader.AddUniform("texsize_x");		glUniform1f(verletShader("texsize_x"),float(texture_size_x));
		verletShader.AddUniform("texsize_y");		glUniform1f(verletShader("texsize_y"),float(texture_size_y));
		verletShader.AddUniform("KsStruct");		glUniform1f(verletShader("KsStruct"),KsStruct);
		verletShader.AddUniform("KdStruct");		glUniform1f(verletShader("KdStruct"),KdStruct);
		verletShader.AddUniform("KsShear");			glUniform1f(verletShader("KsShear"),KsShear);
		verletShader.AddUniform("KdShear");			glUniform1f(verletShader("KdShear"),KdShear);
		verletShader.AddUniform("KsBend");			glUniform1f(verletShader("KsBend"),KsBend);
		verletShader.AddUniform("KdBend");			glUniform1f(verletShader("KdBend"),KdBend);		
		verletShader.AddUniform("inv_cloth_size_x");glUniform1f(verletShader("inv_cloth_size_x"),float(sizeX)/numX);
		verletShader.AddUniform("inv_cloth_size_y");glUniform1f(verletShader("inv_cloth_size_y"),float(sizeY)/numY);		
		verletShader.AddUniform("step");			glUniform2f(verletShader("step"),1.0f/(texture_size_x-1.0f),1.0f/(texture_size_y-1.0f));
	verletShader.UnUse();

	startTime = (float)glutGet(GLUT_ELAPSED_TIME);
	glEnable(GL_DEPTH_TEST);
	int i=0, j=0, count=0;
	int l1=0, l2=0; 
	int v = numY+1;
	int u = numX+1;

	indices.resize( numX*numY*2*3);
 
	X.resize(total_points);
	X_last.resize(total_points); 
	F.resize(total_points);
  
	//fill in positions
	for( j=0;j<=numY;j++) {		 
		for( i=0;i<=numX;i++) {	 
			X[count] = glm::vec4( ((float(i)/(u-1)) *2-1)* hsize, sizeX+1, ((float(j)/(v-1) )* sizeY),1);
			X_last[count] = X[count];
			count++;
		}
	} 

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

	// Setup springs
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
	  
	InitFBO();
	InitVBO(); 
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

void SetOrthographicProjection()
{	
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, 1, 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
}

void ResetPerspectiveProjection() 
{
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

void DrawFullScreenQuad() {
	glBegin(GL_QUADS);
		glVertex2f(0,0);
		glVertex2f(1,0);
		glVertex2f(1,1);
		glVertex2f(0,1);
	glEnd();
}

void RenderGPU() {
	SetOrthographicProjection();
	glViewport(0,0,texture_size_x, texture_size_y);
	CHECK_GL_ERRORS
	for(int i=0;i<NUM_ITER;i++) {
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fboID[writeID]);	
		glDrawBuffers(2, mrt);	
		
		CHECK_GL_ERRORS
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, attachID[2*readID]);
		
		CHECK_GL_ERRORS
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_2D, attachID[2*readID+1]);
		
		glClear(GL_COLOR_BUFFER_BIT);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		verletShader.Use();			 		
			DrawFullScreenQuad();
		verletShader.UnUse();
		//swap read/write pathways
		int tmp = readID;
		readID  = writeID;
		writeID = tmp;
	}
	CHECK_GL_ERRORS
	glFlush();
				
	CHECK_GL_ERRORS						
	ResetPerspectiveProjection();

	//read back the results into the VBO
	glBindFramebuffer(GL_READ_FRAMEBUFFER, fboID[readID]);

	glReadBuffer(GL_COLOR_ATTACHMENT0); 			
	glBindBuffer(GL_PIXEL_PACK_BUFFER, vboID); 			
	glReadPixels(0, 0, texture_size_x, texture_size_y, GL_RGBA, GL_FLOAT, 0); 
						 
	glReadBuffer(GL_NONE); 
	glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);

	//reset default framebuffer
	glBindFramebuffer( GL_FRAMEBUFFER, 0 );
    glReadBuffer(GL_BACK);
    glDrawBuffer(GL_BACK); 
	
	//restore the rendering modes and viewport
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glViewport(0,0,width, height);

	CHECK_GL_ERRORS
 
	renderShader.Use();
		glBindBuffer(GL_ARRAY_BUFFER, vboID);
		glVertexPointer(4, GL_FLOAT, 0,0);
		glEnableClientState(GL_VERTEX_ARRAY);
			//draw plygons
			glUniform4fv(renderShader("color"),1,vWhite);
			glDrawElements(GL_TRIANGLES, (GLsizei)indices.size(), GL_UNSIGNED_SHORT, &(indices[0]));

			//draw points
			glUniform4fv(renderShader("color"),1,vRed);
			glUniform1i(renderShader("selected_index"), selected_index);
			glDrawArrays(GL_POINTS,0, total_points);
			glUniform1i(renderShader("selected_index"), -1);
		glDisableClientState(GL_VERTEX_ARRAY);
	renderShader.UnUse();
 
}
void RenderCPU() {
	 
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
	if ( accumulator >= timeStep ) {	 
        StepPhysics(timeStep );		
        accumulator -= timeStep;
    }

 	//draw polygons
	glColor3f(1,1,1);
	glBegin(GL_TRIANGLES);
	for(i=0;i<indices.size();i+=3) {
		glm::vec3 p1 = vec3(X[indices[i]]);
		glm::vec3 p2 = vec3(X[indices[i+1]]);
		glm::vec3 p3 = vec3(X[indices[i+2]]);
		glVertex3f(p1.x,p1.y,p1.z);
		glVertex3f(p2.x,p2.y,p2.z);
		glVertex3f(p3.x,p3.y,p3.z);
	}
	glEnd();	 

	//draw points	
	glBegin(GL_POINTS);
	for(i=0;i<total_points;i++) {
		glm::vec3 p = vec3(X[i]);
		int is = (i==selected_index);
		glColor3f((float)!is,(float)is,(float)is);
		glVertex3f(p.x,p.y,p.z);
	}
	glEnd();
 
}
void OnRender() {		
	
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

	sprintf_s(info, "FPS: %3.2f, Mode: %s, Frame time (GLUT): %3.4f msecs, Frame time (QP): %3.3f", fps,(current_mode==CPU?"CPU":"GPU"), frameTime, frameTimeQP);
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
	
	switch(current_mode) {
		case CPU:
			RenderCPU();
		break;

		case GPU:
			RenderGPU();
		break;
	} 	
	glutSwapBuffers();
}

void OnShutdown() {
	X.clear();
	X_last.clear();
	F.clear();
	indices.clear();
	springs.clear();

	glDeleteBuffers(1, &vboID);
	glDeleteFramebuffers(2, fboID);
	glDeleteTextures(4, attachID);

	renderShader.DeleteProgram();
	verletShader.DeleteProgram();
}



void IntegrateVerlet(float deltaTime) {
	float deltaTime2 = (deltaTime*deltaTime);
	size_t i=0; 
	

	for(i=0;i<total_points;i++) {		
		glm::vec4 buffer = X[i];
		 
		X[i] = X[i] + (X[i] - X_last[i]) + deltaTime2*glm::vec4(F[i],0);
		  
		X_last[i] = buffer;

		if(X[i].y <0) {
			X[i].y = 0; 
		}
	}
}

inline glm::vec3 GetVerletVelocity(glm::vec3 x_i, glm::vec3 xi_last, float dt ) {
	return  (x_i - xi_last) / dt;
}
void ComputeForces(float dt) {
	size_t i=0;
	 
	for(i=0;i<total_points;i++) {
		F[i] = glm::vec3(0);
		glm::vec3 V = GetVerletVelocity(vec3(X[i]), vec3(X_last[i]), dt);
		//add gravity force
		if(i!=0 && i!=( numX)	)		 
			F[i] += gravity*mass;
		//add force due to damping of velocity
		F[i] += DEFAULT_DAMPING*V;
	}	 

	 
	for(i=0;i<springs.size();i++) {
		glm::vec3 p1 = vec3(X[springs[i].p1]);
		glm::vec3 p1Last = vec3(X_last[springs[i].p1]);
		glm::vec3 p2 = vec3(X[springs[i].p2]);
		glm::vec3 p2Last = vec3(X_last[springs[i].p2]);

		glm::vec3 v1 = GetVerletVelocity(p1, p1Last, dt);
		glm::vec3 v2 = GetVerletVelocity(p2, p2Last, dt);
		
		glm::vec3 deltaP = p1-p2;
		glm::vec3 deltaV = v1-v2;
		float dist = glm::length(deltaP);
		
		float leftTerm = -springs[i].Ks * (dist-springs[i].rest_length);
		float rightTerm = springs[i].Kd * (glm::dot(deltaV, deltaP)/dist);		
		glm::vec3 springForce = (leftTerm + rightTerm)*glm::normalize(deltaP);
			 
		if(springs[i].p1 != 0 && springs[i].p1 != numX)
			F[springs[i].p1] += springForce; 
		if(springs[i].p2 != 0 && springs[i].p2 != numX )
			F[springs[i].p2] -= springForce;
	}
}
 
void ApplyProvotDynamicInverse() {	 
	for(size_t i=0;i<springs.size();i++) { 
		//check the current lengths of all springs
		glm::vec3 p1 = vec3(X[springs[i].p1]);
		glm::vec3 p2 = vec3(X[springs[i].p2]);
		glm::vec3 deltaP = p1-p2;
		float dist = glm::length(deltaP);
		if(dist>springs[i].rest_length) {
			dist -= (springs[i].rest_length);
			dist /= 2.0f;
			deltaP = glm::normalize(deltaP);
			deltaP *= dist;
			if(springs[i].p1==0 || springs[i].p1 ==numX) {
				X[springs[i].p2] += glm::vec4(deltaP,0);
			} else if(springs[i].p2==0 || springs[i].p2 ==numX) {
				X[springs[i].p1] -= glm::vec4(deltaP,0);
			} else { 	
				X[springs[i].p1] -= glm::vec4(deltaP,0);
				X[springs[i].p2] += glm::vec4(deltaP,0);
			}
		}		
	}
}

void OnIdle() {		
	glutPostRedisplay();	
}

void StepPhysics(float dt ) {
	ComputeForces(dt);		
	IntegrateVerlet(dt);   	
}

void OnKey(unsigned char key, int , int) {
	switch(key) {
		case ' ':
			current_mode = (current_mode==CPU)?GPU:CPU;
		break;
	}
	glutPostRedisplay();
}
void main(int argc, char** argv) {
	 
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(width, height);
	glutCreateWindow("GLUT Cloth Demo [Verlet Integration]");

	glutDisplayFunc(OnRender);
	glutReshapeFunc(OnReshape);
	glutIdleFunc(OnIdle);
	
	glutMouseFunc(OnMouseDown);
	glutMotionFunc(OnMouseMove); 
	glutKeyboardFunc(OnKey);
	glutCloseFunc(OnShutdown);

	glewInit();
	InitGL();
	
	glutMainLoop();		
}
