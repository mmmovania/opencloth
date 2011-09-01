//#define USE_SMEM
#define BLOCKSIZE (128 + 2) * (128 + 2)

 
__constant float	KsStruct = 50.75f,KdStruct = -0.25f, 
	KsShear = 50.75f,KdShear = -0.25f,
	KsBend = 50.95f,KdBend = -0.25f;

__constant float3 gravity= (float3)(0.0f,-0.00981f,0.0f); 

int2 getNextNeighbor(int n, float* ks, float* kd) { 
    //structural springs (adjacent neighbors)
    //        o
    //        |
    //     o--m--o
    //        |
    //        o
    if(n<4) {
       *ks = KsStruct;
       *kd = KdStruct;
    }
	if (n == 0)	return (int2)( 1,  0);
	if (n == 1)	return (int2)( 0, -1);
	if (n == 2)	return (int2)(-1,  0);
	if (n == 3)	return (int2)( 0,  1);
	
	//shear springs (diagonal neighbors)
	//     o  o  o
	//      \   /
	//     o  m  o
	//      /   \
	//     o  o  o
	if(n<8) {
       *ks = KsShear;
       *kd = KdShear;
    }
	if (n == 4) return (int2)( 1,  -1);
	if (n == 5) return (int2)( -1, -1);	
	if (n == 6) return (int2)(-1,  1);
	if (n == 7) return (int2)( 1,  1);
	
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
      *ks = KsBend;
      *kd = KdBend;
    }
	if (n == 8)	return (int2)( 2, 0);
	if (n == 9) return (int2)( 0, -2);
	if (n ==10) return (int2)(-2, 0);
	if (n ==11) return (int2)( 0, 2);
}

///////////////////////////////////////////////////////////////////////////////
//! kernel for cloth simulating via verlet integration
//! @param g_odata  memory to process (in and out)
///////////////////////////////////////////////////////////////////////////////
__kernel void verlet(__global float4 * pos_vbo, 
					 __global float4 * g_pos_in, 
					 __global float4 * g_pos_old_in,
					 __global float4 * g_pos_out, 
					 __global float4 * g_pos_old_out, 
					 int texsizeX, 
					 int texsizeY, 
					 float stepX,  
					 float stepY,
					 float damp, 
					 float mass, 
					 float dt, 
					 float inv_cloth_sizeX,
					 float inv_cloth_sizeY)
{
	int2   texsize = (int2)(texsizeX, texsizeY);
	int index = get_global_id(0);
	int ix = index % texsize.x; 
	int iy = index / texsize.x; 

    float ks=0, kd=0;

	float4 pos = g_pos_in[index];    // ensure coalesced read
    float4 pos_old= g_pos_old_in[index];
	float3 vel = (pos.xyz - pos_old.xyz) / dt;
	 
	float2 step = (float2)(stepX, stepY);
	float2 inv_cloth_size = (float2)(inv_cloth_sizeX, inv_cloth_sizeY);
	
	float3 force = gravity*mass + vel*damp;
	  
	if (index==0 || index== (texsize.x - 1.0))
		 mass = 0.0;
	
	for (int k = 0; k < 12; k++)
	{
		int2 coord = getNextNeighbor(k, &ks, &kd);//nextNeigh(k);
		int j = coord.x;
		int i = coord.y;

		if (((iy + i) < 0) || ((iy + i) > (texsize.y - 1)))
			continue;

		if (((ix + j) < 0) || ((ix + j) > (texsize.x - 1)))
			continue;

		int index_neigh = (iy + i) * texsize.x + ix + j;

		float4 pos_neighData = g_pos_in[index_neigh];
		float4 pos_lastData = g_pos_old_in[index_neigh];

		float3 p2 = pos_neighData.xyz;
        float3 p2_last = pos_lastData.xyz;
		float2 coord_neigh = (float2)(ix + j, iy + i) * step;
		
		float rest_length = length((float2)(coord.x*inv_cloth_size.x, coord.y*inv_cloth_size.y));
		 
		 
		float3 v2 = (p2.xyz - p2_last.xyz)/dt;
		float3 deltaP = pos.xyz - p2.xyz;	
		float3 deltaV = vel.xyz - v2.xyz;	 
		float dist = length(deltaP);
				
		float   leftTerm = -ks * (dist-rest_length);
		float  rightTerm = kd * (dot(deltaV, deltaP)/dist);		
		float3 springForce = (leftTerm + rightTerm)*normalize(deltaP);
		force += springForce;	
	}

	float3 acc = (float3)(0, 0, 0);
	if(mass!=0)
		acc = force / mass;

	
	if(pos.y<0)
	   pos.y=0;

	// verlet
	float4 tmp = pos; 
	pos = pos * 2.0 - pos_old + acc * dt * dt;
	pos_old = tmp;
 

	pos_vbo[index] = pos;
	g_pos_out[index] = pos;
	g_pos_old_out[index] = pos_old;
}

