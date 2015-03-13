#OpenCloth
This project was previously hosted on google code. Now that google code will phase out, we have moved the repository to github. 

#Introduction
Introducing an OpenGL based cloth simulation code base. We implement all of the existing cloth simulation algorithms in as simplistic an approach as possible so that learners may know exactly what is needed to get a cloth simulation system up and running in OpenGL with a minimum of fuss. OpenCloth project has been initiated with a view that it may help beginners and researchers alike to implement the basic algorithms for cloth simulation using OpenGL API. It is not intended as another library that you can plugin into your game engine directly. Rather, you can learn from it and then implement a technique or two in you own game/physics engine. With a little bit of effort, it should be straight forward to implement the discussed algorithms on other platforms and Graphics API. Focus is on how to handle the bare minimum required to implement the techniques. Rather than wrapping code into classes, we implement the whole code in a single source file. 

#Publications using OpenCloth
M. Vo, S. G. Narasimhan, and Y. Sheikh, <a target="_blank" href="https://www.cs.cmu.edu/~ILIM/projects/IL/TextIllumSep/papers/CCD14.pdf">Separating Texture and Illumination for Single-Shot Structured Light Reconstruction</a> in The IEEE Conference on Computer Vision and Pattern Recognition (CVPR) Workshops, June 2014. 

Stephen Spinks, <a href="http://www.stephenspinks.com/project.html">Mass Spring Particle Systems</a> Bachelor of Science (Hons.) Final Year Report, University of Glamorgan, May 2012.
 
#Tell us if OpenCloth has helped you in your research
If this project has helped you with your research or project, please let me know. I would like to know more about how the project helped. There are some issues related to the units of simulation used by OpenCloth. Currently, I am trying to consolidate these and hopefully the next release will take care of this.

#Details
Currently. this project contains complete implementations of (in alphabetical order)
  * Co-Rotated Linear FEM
  * Explicit Euler integration
  * Explicit Euler integration with texture mapping and lighting
  * Explicit Euler integration with wind
  * Implicit Explicit (IMEX) method
  * Implicit integration (Baraff & Witkin's model)
  * Implicit Euler integration
  * Meshless FEM 
  * Position based dynamics
  * Semi-Implicit integration (Symplectic Euler)
  * Verlet integration
  * Verlet integration on CUDA, GLSL (using GPGPU technique) and OpenCL
  * WebGL port of Explicit Euler Integration 

#Live WebGL Demos
<a target="_blank" href="http://opencloth.googlecode.com/svn/trunk/OpenCloth_WebGL/WebGLOpenCloth.html">Live WebGL Demo 1</a>

<a target="_blank" href="http://opencloth.googlecode.com/svn/trunk/OpenCloth_WebGL/WebGLOpenClothTextured.html">Live WebGL Demo 2</a>

#Demo Videos
This section will add more recent videos that reflect the recent contents of OpenCloth.

#Explicit Euler integration
[![OpenCloth Demo 1](http://img.youtube.com/vi/5MuzlGmLngY/0.jpg)](http://www.youtube.com/watch?v=5MuzlGmLngY)

#Manipulation
[![OpenCloth Demo 2](http://img.youtube.com/vi/2E7h38U5-as/0.jpg)](http://www.youtube.com/watch?v=2E7h38U5-as) 
