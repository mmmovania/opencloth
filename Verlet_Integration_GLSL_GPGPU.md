# Introduction #
We may perform the Verlet integration on the GPU. There are three approaches to do this.
  1. Using the vertex shader
  1. Using the GPGPU based approach which exploits the fragment shader.
  1. Using my approach (which I will describe once my paper is published :P)
I have implemented the second approach and my reference was this wonderful paper by Joachim et. al http://wwwcg.in.tum.de/Research/data/Publications/simpra05.pdf
. In case any of the things are unclear, do refer to the original article for details.

# Details #
If you read through the cited article, you would see two approaches mentioned in it,
  1. Point centric approach and
  1. Edge centric approach
I went with the point-centric approach for my implementation.
## Data structure ##
Since this approach uses the fragment shader, it uses the standard GPGPU approach whereby a full screen quad is rendered and then a fragment shader is invoked. Force calculation and constraints are applied per pixel and the results are written to an offscreen texture. To aid in this, I use the ping pong approach, which involves writing to a render target (draw framebuffer) texture while reading from another render target (read framebuffer). First the FBO and vertex buffers are setup to contain the cloth vertex positions. For Verlet integration, two buffers are needed so that the previous position may also be stored. Initially the current positions are written to the previous position texture. In the render function, the quad is rendered and then the fragment shader is invoked. To give an idea of whats going on in the fragment shader, have a look at the following pseudocode.
```
Texture positionTexture;       //current positions
Texture prevPositionTexture;   //previous positions
void main() {
   Get current position by looking up the position texture
   Get previous position by looking up the prev. position texture
   Calculate initial velocity using Verlet integration (v = (x_i-x_last)/dt)
   Calculate index from the uv coordinates and texture size
   Initialize force with gravity and damping the initial velocity   
   for each neighbor n
      Find coordinates of n
      Get position of n and compute its spring force 
      Add spring force to the total force 
   end for

   Calculate acceleration using a=F/m;
   Calculate new position using Verlet integration (x_i=2x_i-x_last + a*t*t)
   Apply collision constraints
   Write positions to shader outputs
}
```
That's it. The only tricky part in this code is the actual setup of FBO and the ping pong strategy. The rest of the things are all standard Verlet implementation. I hope this document clarifies the GLSL based code for Verlet integration.