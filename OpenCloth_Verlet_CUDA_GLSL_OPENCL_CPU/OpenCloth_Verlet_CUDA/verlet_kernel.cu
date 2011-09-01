#include <cutil_math.h>
#include <math_constants.h>

//#define USE_SMEM
#define BLOCKSIZE (128 + 2) * (128 + 2)

 
const __device__ float	KsStruct = 50.75f,KdStruct = -0.25f, 
	KsShear = 50.75f,KdShear = -0.25f,
	KsBend = 50.95f,KdBend = -0.25f;



__device__ int2 getNextNeighbor(int n, float& ks, float& kd) { 
    //structural springs (adjacent neighbors)
    //        o
    //        |
    //     o--m--o
    //        |
    //        o
    if(n<4) {
       ks = KsStruct;
       kd = KdStruct;
    }
	if (n == 0)	return make_int2( 1,  0);
	if (n == 1)	return make_int2( 0, -1);
	if (n == 2)	return make_int2(-1,  0);
	if (n == 3)	return make_int2( 0,  1);
	
	//shear springs (diagonal neighbors)
	//     o  o  o
	//      \   /
	//     o  m  o
	//      /   \
	//     o  o  o
	if(n<8) {
       ks = KsShear;
       kd = KdShear;
    }
	if (n == 4) return make_int2( 1,  -1);
	if (n == 5) return make_int2( -1, -1);	
	if (n == 6) return make_int2(-1,  1);
	if (n == 7) return make_int2( 1,  1);
	
	//bend spring (adjacent neighbors 1 node away)
	//
	//o   o   o   o   o
	//        | 
	//o   o   |   o   o
	//        |   
	//o-------m-------o
	//        |  
	//o   o   |   o   o
	//        |
	//o   o   o   o   o 
	if(n<12) {
       ks = KsBend;
       kd = KdBend;
    }
	if (n == 8)	return make_int2( 2, 0);
	if (n == 9) return make_int2( 0, -2);
	if (n ==10) return make_int2(-2, 0);
	if (n ==11) return make_int2( 0, 2);
}

///////////////////////////////////////////////////////////////////////////////
//! kernel for cloth simulating via verlet integration
//! @param g_odata  memory to process (in and out)
///////////////////////////////////////////////////////////////////////////////
__global__ void verlet(	float4 * pos_vbo, float4 * g_pos_in, float4 * g_pos_old_in, float4 * g_pos_out, float4 * g_pos_old_out, 
							int2 texsize, float2 step,  float damp, float mass, float dt, float2 inv_cloth_size)
{
    uint index = blockIdx.x * blockDim.x + threadIdx.x;
	int ix = index % texsize.x; 
	int iy = index / texsize.x; 

	//g_pos[index] = make_float4(threadIdx.x, blockIdx.x, blockDim.x, index);
	//return;
    float ks=0, kd=0;
#ifdef USE_SMEM
    __shared__ float4 smem_pos[BLOCKSIZE];
    __shared__ float4 smem_last_pos[BLOCKSIZE];

	int ix_smem = threadIdx.x % texsize.x;  
	int iy_smem = threadIdx.x / texsize.x; 

	smem_pos[threadIdx.x] = g_pos_in[index]; 
	smem_last_pos[threadIdx.x] = g_pos_old_in[index]; 

	for (int k = 0; k < 12; k++)
	{
		int2 coord = getNextNeighbor(k, ks, kd);
		int j = coord.x;
		int i = coord.y;

		if (((iy_smem + i) < 0) || ((iy_smem + i) > (texsize.x - 1)))
			continue;

		if (((ix_smem + j) < 0) || ((ix_smem + j) > (texsize.x - 1)))
			continue;

		int index_neigh_smem = (iy_smem + i) * texsize.x + ix_smem + j;
		int index_neigh = (iy + i) * texsize.x + ix + j;

		smem_pos[index_neigh_smem] = g_pos_in[index_neigh]; 
		smem_last_pos[index_neigh_smem] = g_pos_old_in[index_neigh]; 
	}

	__syncthreads();

	volatile float4 posData = smem_pos[threadIdx.x];    // ensure coalesced read
    volatile float4 posOldData = smem_last_pos[threadIdx.x];
#else
	volatile float4 posData = g_pos_in[index];    // ensure coalesced read
    volatile float4 posOldData = g_pos_old_in[index];
#endif


    float3 pos = make_float3(posData.x, posData.y, posData.z);
    float3 pos_old = make_float3(posOldData.x, posOldData.y, posOldData.z);
	float3 vel = (pos - pos_old) / dt;
	 
	const float3 gravity=make_float3(0.0f,-0.00981f,0.0f); 
	float3 force = gravity*mass + vel*damp;
	  
	if (index==0 || index== (texsize.x - 1.0))
		 mass = 0.0;

	
	for (int k = 0; k < 12; k++)
	{
		int2 coord = getNextNeighbor(k, ks, kd);//nextNeigh(k);
		int j = coord.x;
		int i = coord.y;

#ifdef USE_SMEM
		if (((iy_smem + i) < 0) || ((iy_smem + i) > (texsize.x - 1)))
			continue;

		if (((ix_smem + j) < 0) || ((ix_smem + j) > (texsize.x - 1)))
			continue;

		int index_neigh_smem = (iy_smem + i) * texsize.x + ix_smem + j;

		volatile float4 pos_neighData = smem_pos[index_neigh_smem];
		volatile float4 pos_lastData = smem_last_pos[index_neigh_smem];
		
#else
		if (((iy + i) < 0) || ((iy + i) > (texsize.y - 1)))
			continue;

		if (((ix + j) < 0) || ((ix + j) > (texsize.x - 1)))
			continue;

		int index_neigh = (iy + i) * texsize.x + ix + j;

		volatile float4 pos_neighData = g_pos_in[index_neigh];
		volatile float4 pos_lastData = g_pos_old_in[index_neigh];
#endif
		float3 p2 = make_float3(pos_neighData.x, pos_neighData.y, pos_neighData.z);
        float3 p2_last = make_float3(pos_lastData.x, pos_lastData.y, pos_lastData.z);
		float2 coord_neigh = make_float2(ix + j, iy + i) * step;
		
		float rest_length = length(make_float2(coord.x*inv_cloth_size.x, coord.y*inv_cloth_size.y));
		 
		 
		float3 v2 = (p2- p2_last)/dt;
		float3 deltaP = pos - p2;	
		float3 deltaV = vel - v2;	 
		float dist = length(deltaP);
				
		float   leftTerm = -ks * (dist-rest_length);
		float  rightTerm = kd * (dot(deltaV, deltaP)/dist);		
		float3 springForce = (leftTerm + rightTerm)*normalize(deltaP);
		force += springForce;	
	}

	float3 acc = make_float3(0, 0, 0);
	if(mass!=0)
		acc = force / mass;

	
	if(pos.y<0)
	   pos.y=0;

	// verlet
	float3 tmp = pos; 
	pos = pos * 2 - pos_old + acc * dt * dt;
	pos_old = tmp;

	syncthreads();

	pos_vbo[index] = make_float4(pos.x, pos.y, pos.z, posData.w);
	g_pos_out[index] = make_float4(pos.x, pos.y, pos.z, posData.w);
	g_pos_old_out[index] = make_float4(pos_old.x, pos_old.y, pos_old.z, posOldData.w);

}

