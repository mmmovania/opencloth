#include <cutil_inline.h>
#include <cutil_math.h>

#include <cassert>

float4* X[2];
float4* X_last[2];

float4 * X_in, *X_out;
float4 * X_last_in, *X_last_out;

extern int readID, writeID;
__global__ void verlet(	float4 * pos_vbo, float4 * g_pos_in, float4 * g_pos_old_in, float4 * g_pos_out, float4 * g_pos_old_out, 
							int2 texsize, float2 step,  float damp, float mass, float dt, float2 inv_cloth_size);

void InitCUDA(const unsigned int size) {
	const unsigned int num_threads = size;
	const unsigned int mem_size = sizeof(float4) * num_threads;
	
	// allocate device memory for float4 version
	cutilSafeCall(cudaMalloc((void**) &X[0], mem_size));	// positions
	cutilSafeCall(cudaMalloc((void**) &X[1], mem_size));	// positions
	cutilSafeCall(cudaMalloc((void**) &X_last[0], mem_size));	// old positions
	cutilSafeCall(cudaMalloc((void**) &X_last[1], mem_size));	// old positions		
}

void ShutdownCUDA()
{
	// cleanup memory
	if (X[0] != NULL) 
	{
		cutilSafeCall(cudaFree(X[0]));
		cutilSafeCall(cudaFree(X[1]));
		X[0] = NULL;
		X[1] = NULL;
	}

	if (X_last[0] != NULL)
	{
		cutilSafeCall(cudaFree(X_last[0]));
		cutilSafeCall(cudaFree(X_last[1]));
		X_last[0] = NULL;
		X_last[1] = NULL;
	}
}

void computeGridSize(uint n, uint blockSize, uint &numBlocks, uint &numThreads)
{
    numThreads = min(blockSize, n); 
    numBlocks = (n % numThreads != 0) ? (n / numThreads + 1) : (n / numThreads);
}  

void UploadCUDA(float * positions, float * positions_old, const int size)
{
	static bool start = true;

	assert(X[0] != NULL); 
	assert(X_last[0] != NULL); 

	const unsigned int num_threads = size;
	const unsigned int mem_size = sizeof(float4) * num_threads;

	X_in  = X[readID];	
	X_out = X[writeID];
	X_last_in  = X_last[readID];	
	X_last_out = X_last[writeID];
	
	if (start)
	{
		cutilSafeCall(cudaMemcpy(X_in, positions,  mem_size, cudaMemcpyHostToDevice));
		cutilSafeCall(cudaMemcpy(X_last_in, positions_old, mem_size, cudaMemcpyHostToDevice));
		cutilCheckMsg("Cuda memory copy host to device failed.");
		start=false;
	} 

	int tmp=readID;
	readID = writeID;
	writeID=tmp;
}

void VerletCUDA(float4 * pos_vbo, int2 texsize, float2 step, const float & damp, const float & mass, float dt, float2 inv_cloth_size)
{   
	// setup execution parameters 
	uint numThreads, numBlocks;
	uint numParticles = texsize.x*texsize.y;

	computeGridSize(numParticles, 256, numBlocks, numThreads);

//	printf("%3d particles, %3d blocks, %3d threads\n", numParticles, numBlocks, numThreads);

	// execute the kernel
	//	printf("numParticles: %d,   numThreads: %d   numBlocks: %d\n", numParticles, numThreads, numBlocks);
	verlet<<< numBlocks, numThreads >>>(pos_vbo, X_in, X_last_in, X_out, X_last_out, texsize, step, damp, mass, dt, inv_cloth_size);

	// stop the CPU until the kernel has been executed
	cudaThreadSynchronize();

	// check if kernel execution generated and error
	cutilCheckMsg("Cuda kernel execution failed.");
}
