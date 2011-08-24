
uniform sampler2D X;		//current position
uniform sampler2D X_last;	//previous position
uniform float dt;			//timeStep
uniform float DEFAULT_DAMPING;
uniform float inv_mass;
uniform float texsize_x;
uniform float texsize_y;
uniform float KsStruct,KdStruct,KsShear,KdShear,KsBend,KdBend;
uniform float inv_cloth_size_x;
uniform float inv_cloth_size_y;
uniform vec2  step;

vec3 gravity = vec3(0, -9.81, 0);

vec2 getNextNeighbor(int n, out float ks, out float kd) { 
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
	if (n == 0)	return vec2( 1,  0);
	if (n == 1)	return vec2( 0, -1);
	if (n == 2)	return vec2(-1,  0);
	if (n == 3)	return vec2( 0,  1);
	
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
	if (n == 4) return vec2( 1,  -1);
	if (n == 5) return vec2( -1, -1);	
	if (n == 6) return vec2(-1,  1);
	if (n == 7) return vec2( 1,  1);
	
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
	if (n == 8)	return vec2( 2, 0);
	if (n == 9) return vec2( 0, -2);
	if (n ==10) return vec2(-2, 0);
	if (n ==11) return vec2( 0, 2);
}

 

void main() {
	
	vec3 x_i	= texture2D(X, gl_TexCoord[0].st).xyz;
	vec3 x_last = texture2D(X_last, gl_TexCoord[0].st).xyz;
	vec3 vel	= (x_i - x_last) / dt;	// calc. velocity according to verlet integration
	
	//vec2 step = vec2(1.0,1.0) / (texsize_x - 1.0, texsize_y-1.0);
	
	float ix = floor(gl_TexCoord[0].s * texsize_x);
	float iy = floor(gl_TexCoord[0].t * texsize_y);
	float index = iy * texsize_x + ix;
		
	float inverse_mass = inv_mass;
	float ks=  0.0, kd= 0.0;
	
	 
	if (index==0 || index== (texsize_x - 1.0))
		 inverse_mass = 0.0;
 
	vec3 force = gravity*inverse_mass + vel*DEFAULT_DAMPING;
 
	for (int k = 0; k < 12; k++)
	{ 	    
		vec2 coord = getNextNeighbor(k, ks, kd);
		float j = coord.x;
		float i = coord.y;

		if (((iy + i) < 0.0) || ((iy + i) > (texsize_y - 1.0)))
			continue;

		if (((ix + j) < 0.0) || ((ix + j) > (texsize_x - 1.0)))
			continue;

		vec2 coord_neigh = vec2(ix + j, iy + i) * step;
		
		float rest_length = length(coord)*inv_cloth_size_x;
		//float rest_length =  max(abs(coord.x*inv_cloth_size_x), abs(coord.y*inv_cloth_size_y));

		vec3 p2 = texture2D(X, coord_neigh).xyz;
		 
		vec3 deltaP = x_i - p2;		 
		float dist = length(deltaP);
		
		float   leftTerm = -ks * (dist-rest_length);
		float  rightTerm = kd * (dot(vel, deltaP)/dist);		
		vec3 springForce = (leftTerm + rightTerm)*normalize(deltaP);
		force += springForce;		
	} 
	
	
	vec3 acc;
	acc = force*inverse_mass;

	if(x_i.y<0)
	   x_i.y=0;
	   
	// verlet
	vec3 tmp = x_i;
	x_i = x_i * 2.0 - x_last + acc * dt * dt;
	x_last = tmp;
		
	gl_FragData[0] = vec4(x_i,1.0);
	gl_FragData[1] = vec4(x_last,0.0);	
}