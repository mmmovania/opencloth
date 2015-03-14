# Introduction #

Introducing an OpenGL based cloth simulation code base. We implement all of the existing cloth simulation algorithms in as simplistic an approach as possible so that learners may know exactly what is needed to get a cloth simulation system up and running in OpenGL with a minimum of fuss. OpenCloth project has been initiated with a view that it may help beginners and researchers alike to implement the basic algothims for cloth simulation using OpenGL API. It is not intended as another library that you can plugin into your game engine directly. Rather, you can learn from it and then implement a technique or two in you own game/physics engine. With a little bit of effort, it should be straight forward to implement the discussed algorithms on other platforms and Graphics API.

While at the current initial stage, we are just releasing a bunch of source codes, later on, we will try to implement these in the form of a platform independent library.

# Details #

Current codes contain the complete implementations of (in alphabetical order)
  1. Co-Rotated Linear FEM
  1. Explicit Euler integration
  1. Explicit Euler integration with texture mapping and lighting
  1. Explicit Euler integration with wind
  1. Implicit Explicit (IMEX) method
  1. Implicit integration (Baraff & Witkin's model)
  1. Implicit Euler integration
  1. Meshless FEM
  1. Position based dynamics
  1. Semi-Implicit integration (Symplectic Euler)
  1. Verlet integration
  1. Verlet integration on CUDA, GLSL (using GPGPU technique) and OpenCL

## What is not offered in the current version ##
Currently, the collision detection/response is only with an arbitrarily oriented ellipsoid and the naive collision to the floor. More primitives like spheres, cubes, convex hulls, arbitrary meshes are in the pipeline and will be added as time permits. Spatial indexing datastructures and Kd-Trees for fast neighbor searches are also in the pipeline.

## External dependencies ##
The codes depend on the freeglut, glew and glm libraries. For code using CUDA/OpenCL, the latest version of the GPU computing sdk is needed. These are not included with this project and you would need to download them separately.

## Understanding the code ##
If the code does not make sense, I would suggest a read of the SIGGRAPH course notes (Realtime Physics: www.matthiasmueller.info/realtimephysics/coursenotes.pdf ).

For bug reports/suggestions/queries email to me at mova0002@e.ntu.edu.sg.

## To Do List ##
Some features that I will be adding soon.
  * Hierarchical position based dynamics
  * Collision detection and response (spheres,cubes,convex hulls etc.)
  * Demo application encapsulating all algorithms
  * Extend to OpenGL3.3 core profile
  * Platform independent library (C++ class hierarchy)
  * Soft bodies (Tetrahedra/Hexahedra FEM)