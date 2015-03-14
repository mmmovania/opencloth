# Introduction #
The initial IMEX implementation in OpenCloth was fundametally flawed thanks to Saggita for spotting this to me. Honestly, I simply copied content from the algorithm given in the Reference [1](1.md). To reimplement it correctly, I had to dive deep into the paper to understand what it is trying to do. So here is my understanding about this awesome method.

# Details #
We know that there are two basic classes of integrations for numercial simulations namely implicit and explicit. The former are more stable but require calculation of a large matrix that is solved using PCG (Preconditioned Conjugate Gradient or similar methods. The latter are simple but unstable. IMEX are a class of integration that give a nice middle ground between implicit and explicit methods.

# IMEX Method #
I will now detail the IMEX method for cloth simulation as detailed in Reference [1](1.md). The first thing that we need is the Jacobian matrix (J) which stores the partial derivatives of the forces at each mass point. The reference [1](1.md) writes it as Hessian matrix (H). This is a very large matrix N\*N for a cloth mesh with N vertices.
In case of 3D, [1](1.md) defines H as
```
if(i != j)
   Hij=kij
else
   Hij=-sum_i kij
```
This matrix is easily obtained by using the spring stiffness values (kij).

# Calculating the W matrix using Hessian matrix (H) #
The core equation is given on page 3 Eq (3) which is as follows
```
V^n+1 - V^n = (I - (dt^2)/m * H)^-1 (F^n + dt*H*v^n)*dt/m
```
The part `(I - (dt^2)/m * H)` is combined to get a matrix W. The inverse of this matrix is used as an approximation of implicit integration during the predictor stage as detailed below. Since W is an N\*N matrix (N is the total number of mass points), we find the inverse using the Armadillo C++ library which gives fast and convenient matrix routines. Like the original author, I do not store 3N\*3N values rather just the N\*N values to save space and processing time.

# Force Calculation and IMEX Integration #
Reference [1](1.md) splits the forces into two parts: a linear part and a non-linear part. The linear part uses the Hessian matrix defined earlier to approximate implicit integration. This integration provides what is called a predictor for velocity and position. The non-linear part is solved by using a simply correction method that tries to balance the rotation error caused by the non-linear forces. This is implemented in the correcter code.

# Internal Force Calculation #
The original paper by Desbrun et al. [1](1.md) gives the following formula for the calculation of internal force.
```
vec3 F_linear = -springs[i].Ks * (dist-springs[i].rest_length) * 
                normalize(deltaP);
vec3 F_nonlinear = -springs[i].Ks * deltaV * timeStep ;
vec3 springForce = (F_linear + F_nonlinear ) ;
```
This unfortunately gives a severly damped result. I modified the spring force to following which gives a much better result
```
float leftTerm  = -springs[i].Ks * (dist-springs[i].rest_length);
float rightTerm = -springs[i].Kd * (dot(deltaV, deltaP)/dist);		
vec3 springForce = (leftTerm + rightTerm)*normalize(deltaP);
```

# Implementation #
OK now the fun part. How to put all this in code. I will just highlight the crux. The rest of the things should be easily understood by following the code. We first precompute W as follows:
```
//calculate off-diagonal and diagonal terms of H matrix (see Eq. 4 on page 4)
for(size_t j=0;j<total_points;++j) {
   float sum = 0;
   for(size_t i=0;i<total_points;++i) {
      if(i!=j) {
         H[i+j*total_points] = KsStruct;
	 sum += KsStruct;
      }
      H[j+j*total_points] = sum;
   }
}

//Calculate the inverse using Armadillo library
arma::mat T(total_points, total_points); 
for(int j=0;j<total_points;++j) {
   for(int i=0;i<total_points;++i) {
      int index = i+j*total_points; 
      T(j,i) = (1 - dt2_m*H[index]);
   }
}
T = T.i();
	
for(int j=0;j<total_points;++j) {
   for(int i=0;i<total_points;++i) {
      int index = i+j*total_points;
      W[index] = T(j,i);
   }
}
```
Here is the code for predictor
```
glm::vec3 Xg = glm::vec3(0);
glm::vec3 delTau = glm::vec3(0);

//Predictor
for(i=0;i<total_points;i++) {
   glm::vec3 sum=glm::vec3(0);
   glm::vec3 F_filtered=glm::vec3(0);
   for(j=0;j<total_points;++j)
      if(j!=i)
         sum  += F[j]*W[i*total_points + j];
   F_filtered = sum;
   delTau += glm::cross(F_filtered, X[i]);
   Xg += X[i];

   V[i] += ((F[i] + F_filtered)*deltaTimeMass);				
		 
   X_new[i] = X[i] + deltaTime*V[i];		
}
Xg /= total_points; 
```

The corrector code is relatively straightforward. It simply uses the corrected force to estimate the new position as shown below
```
//Corrector
glm::vec3 F_corrected = glm::vec3(0);
for(i=0;i<total_points;i++) {
   if(i!=0 && i!=( numX) ) {
      F_corrected = glm::cross( (Xg - X[i]), delTau) * deltaTime;
      X_new[i] += F_corrected*deltaT2Mass; 
   }
}
```
After the calculation of the new corrected position, we update the velocities and positions
```
void UpdateRealVelocityAndPosition() {
   static float invDeltaTime = 1.0f/timeStep;
   for(size_t i=0;i<total_points;i++) { 
      V[i] = (X_new[i]-X[i])*invDeltaTime;	
      V[i] *= GLOBAL_DAMPING;
      X[i] = X_new[i];	
      if(X[i].y <0) {
         X[i].y = 0; 
      }
   }
}
```

So to sum up, the physics update proceeds as follows.
```
void StepPhysics(float dt ) {
   ComputeForces();		
   IntegrateSemiImplicit(dt);
   EllipsoidCollision();
   for(int i=0;i<NUM_ITER;++i)
      ApplyProvotDynamicInverse();	
   UpdateRealVelocityAndPosition();
}
```

# Reference #
[1](1.md) Interactive animation of structured deformable objects by Desbrun et al. url: http://www.multires.caltech.edu/pubs/GI99.pdf