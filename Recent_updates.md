# Introduction #

This page is higlighting the new changes that will be done in the code base.


# Details #
  * A slight update, I am currently working on removing the issues, I have already removed the IMEX integration problem. I will check in the updated code soon. Now I will move on to issue related to Implicit Integration.
  * Meanwhile I have also worked on self collision with the basic point triangle and edge-edge collision test done.
  * I will be doing a basic self collision demo soon using the Provot seminal paper (http://graphics.stanford.edu/courses/cs468-02-winter/Papers/Collisions_vetements.pdf) as a start.
  * Since collision detection requires a lot of neighbor search, I will be adding support for collision structures like spatial hashing and KdTrees to reduce the neighbor search.
  * Once these are ready, we will see how to add suport for more current state-of-the-art self collision methods for cloth self collision.