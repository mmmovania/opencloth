//OpenCL variables
#include <oclUtils.h>
cl_platform_id cpPlatform;
cl_context cxGPUContext;
cl_command_queue cqCommandQueue;
cl_program	hProgram;
cl_kernel	hKernel;
 
cl_mem		hVbo;


cl_mem X[2];
cl_mem X_last[2];
 

extern int readID, writeID;


void createProgramAndKernel(const char * clSourcefile, const char * clKernelName, cl_program & cpProgram, cl_kernel & ckKernel)
{
    cl_int ciErrNum;

	// Program Setup
    size_t program_length;
    char * source = oclLoadProgSource(clSourcefile, "", &program_length);
    oclCheckError(source != NULL, shrTRUE);

    // create the program
    cpProgram = clCreateProgramWithSource(cxGPUContext, 1,(const char **) &source, &program_length, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);
    free(source);

    // build the program
    ciErrNum = clBuildProgram(cpProgram, 0, NULL, NULL/*"-cl-fast-relaxed-math"*/, NULL, NULL);
    if (ciErrNum != CL_SUCCESS)
    {
        // write out standard error, Build Log and PTX, then cleanup and exit
        oclLogBuildInfo(cpProgram, oclGetFirstDev(cxGPUContext));
        oclLogPtx(cpProgram, oclGetFirstDev(cxGPUContext), "oclVerlet.ptx");
//      Cleanup(EXIT_FAILURE); 
		shrLogEx(LOGBOTH | CLOSELOG, 0, "GPGPU cloth exiting...\nPress <Enter> to Quit\n");
        shrLogEx(LOGBOTH | ERRORMSG, ciErrNum, STDERROR);
		exit(EXIT_FAILURE);
    }

    // create the kernel
    ckKernel = clCreateKernel(cpProgram, clKernelName, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);
}

void CreateCL_VBO(const int size) {
	cl_int ciErrNum;
	hVbo=clCreateBuffer(cxGPUContext, CL_MEM_WRITE_ONLY, size, NULL, &ciErrNum);    
    oclCheckError(ciErrNum, CL_SUCCESS);
}

void InitOpenCL(const unsigned int size, int texture_size_x, int texture_size_y,float stepX, float stepY, float damp, float mass, float dt, float inv_cloth_sizeX, float inv_cloth_sizeY) {
	//Setup OpenCL context
	//get device ID
	cl_device_id cdDevice;
	cl_int ciErrNum;

    oclCheckError(oclGetPlatformID(&cpPlatform), CL_SUCCESS);
	oclCheckError(clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 1, &cdDevice, NULL), CL_SUCCESS);
	cl_context_properties props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)cpPlatform, 0};

	cxGPUContext = clCreateContext(0, 1, &cdDevice, NULL, NULL, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);
	
	oclPrintDevInfo(LOGBOTH, cdDevice);
	
	cqCommandQueue = clCreateCommandQueue(cxGPUContext, cdDevice, 0, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);

	oclPrintDevName(LOGBOTH, cdDevice);

	const unsigned int num_threads = size;
	const unsigned int mem_size = 4*sizeof(float) * num_threads;
	
	//create 4 buffers for current and previous positions
	X[0] = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE, mem_size , NULL, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);
	X_last[0] = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE, mem_size , NULL, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);

	X[1] = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE, mem_size , NULL, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);	
	X_last[1] = clCreateBuffer(cxGPUContext, CL_MEM_READ_WRITE, mem_size , NULL, &ciErrNum);
    oclCheckError(ciErrNum, CL_SUCCESS);	

	const char* clSourcefile = "verlet_kernel.cl";
	const char* clKernelName = "verlet";
	createProgramAndKernel(clSourcefile, clKernelName, hProgram, hKernel);

	CreateCL_VBO(mem_size);

	ciErrNum |= clSetKernelArg(hKernel, 0, sizeof(cl_mem), (void *)&(hVbo));
	
	ciErrNum |= clSetKernelArg(hKernel, 5, sizeof(int), (void *)&(texture_size_x));
	ciErrNum |= clSetKernelArg(hKernel, 6, sizeof(int), (void *)&(texture_size_y));
	ciErrNum |= clSetKernelArg(hKernel, 7, sizeof(float), (void *)&(stepX));
	ciErrNum |= clSetKernelArg(hKernel, 8, sizeof(float), (void *)&(stepY));
	ciErrNum |= clSetKernelArg(hKernel, 9, sizeof(float), (void *)&(damp));
	ciErrNum |= clSetKernelArg(hKernel, 10, sizeof(float), (void *)&(mass));
	ciErrNum |= clSetKernelArg(hKernel, 11, sizeof(float), (void *)&(dt));
	ciErrNum |= clSetKernelArg(hKernel, 12, sizeof(float), (void *)&(inv_cloth_sizeX));
	ciErrNum |= clSetKernelArg(hKernel, 13, sizeof(float), (void *)&(inv_cloth_sizeY));

	oclCheckError(ciErrNum, CL_SUCCESS);


}

void ShutdownOpenCL()
{
	cl_int ciErrNum = CL_SUCCESS;
   
	// cleanup memory
	if (X[0] != NULL) 
	{ 
		ciErrNum |= clReleaseMemObject(X[0]);
		ciErrNum |= clReleaseMemObject(X[1]);
		X[0] = NULL;
		X[1] = NULL;
	}

	if (X_last[0] != NULL)
	{ 
		ciErrNum |= clReleaseMemObject(X_last[0]);
		ciErrNum |= clReleaseMemObject(X_last[1]);

		X_last[0] = NULL;
		X_last[1] = NULL;
	}
	oclCheckError(ciErrNum, CL_SUCCESS);

	ciErrNum |= clReleaseMemObject(hVbo);
	oclCheckError(ciErrNum, CL_SUCCESS);

	ciErrNum |= clReleaseProgram(hProgram);
	oclCheckError(ciErrNum, CL_SUCCESS);

	ciErrNum |= clReleaseKernel(hKernel);
	oclCheckError(ciErrNum, CL_SUCCESS);

	ciErrNum |= clReleaseCommandQueue(cqCommandQueue);
	oclCheckError(ciErrNum, CL_SUCCESS);

    ciErrNum |= clReleaseContext(cxGPUContext);
    oclCheckError(ciErrNum, CL_SUCCESS);
}

  

void UploadOpenCL(float * positions, float * positions_old, const int size)
{
	static bool start = true;
	const unsigned int mem_size = 4*sizeof(float) * size;
	
	cl_int ciErrNum=CL_SUCCESS;
	
	assert(X[readID] != NULL); 
	assert(X_last[readID] != NULL); 

	if (start)
	{
		//Copy the current positions
		ciErrNum |= clEnqueueWriteBuffer(cqCommandQueue, X[readID], CL_TRUE, 0, mem_size, positions, 0, NULL, NULL);
		oclCheckError(ciErrNum, CL_SUCCESS);
		ciErrNum |= clEnqueueWriteBuffer(cqCommandQueue, X_last[readID], CL_TRUE, 0, mem_size, positions_old, 0, NULL, NULL);
		oclCheckError(ciErrNum, CL_SUCCESS); 
		start=false;
	}

	ciErrNum |= clSetKernelArg(hKernel, 1, sizeof(cl_mem), (void *)&X[readID]);
	ciErrNum |= clSetKernelArg(hKernel, 2, sizeof(cl_mem), (void *)&X_last[readID]);
	ciErrNum |= clSetKernelArg(hKernel, 3, sizeof(cl_mem), (void *)&X[writeID]);
	ciErrNum |= clSetKernelArg(hKernel, 4, sizeof(cl_mem), (void *)&X_last[writeID]);
	
	oclCheckError(ciErrNum, CL_SUCCESS);

	int tmp=readID;
	readID = writeID;
	writeID=tmp;
}

size_t uSnap(size_t a, size_t b)
{
	return ((a % b) == 0) ? a : (a - (a % b) + b);
}

void VerletOpenCL(int texsizeX, int texsizeY)
{   
	// setup execution parameters 
	//uint numThreads, numBlocks;
	int numParticles = texsizeX*texsizeY;
	
	//computeGridSize(numParticles, 256, numBlocks, numThreads);

	size_t _szLocalWorkSize = 256;
	size_t _szGlobalWorkSize = uSnap(numParticles, _szLocalWorkSize);

//	printf("%3d particles, %3d blocks, %3d threads\n", numParticles, numBlocks, numThreads);

	// execute the kernel
	//	printf("numParticles: %d,   numThreads: %d   numBlocks: %d\n", numParticles, numThreads, numBlocks);
	//verlet<<< numBlocks, numThreads >>>(pos_vbo, X_in, X_last_in, X_out, X_last_out, texsize, step, damp, mass, dt, inv_cloth_size);
	cl_int ciErrNum = clEnqueueNDRangeKernel(cqCommandQueue, hKernel, 1, NULL, &_szGlobalWorkSize, &_szLocalWorkSize, 0,0,0 );
    oclCheckError(ciErrNum, CL_SUCCESS);

	// stop the CPU until the kernel has been executed
	//cudaThreadSynchronize();
	ciErrNum |= clFinish(cqCommandQueue);

	// check if kernel execution generated and error
	oclCheckError(ciErrNum, CL_SUCCESS);
}

void ReadBuffer(float* ptr, int size) {
	cl_int ciErrNum = clEnqueueReadBuffer(cqCommandQueue, hVbo, CL_TRUE, 0, size, ptr, 0, NULL, NULL);	
	oclCheckError(ciErrNum, CL_SUCCESS);	
}