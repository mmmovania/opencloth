# Introduction #

This document will illustrate the steps needed to add external forces like Wind and collision detection to the cloth simulation.


# Explanation #
## Adding wind force ##
Adding wind should be pretty straightforward. Just add another force (F\_wind) in the external force calculation function. Currently, it is something like this
```
void computeForce() {
   for i=0 to total_positions
      F[i] = gravity + Kd*V[i];
} 
```
you may add in F\_wind like this
```
void computeForce() {
   for i=0 to total_positions
      F[i] = gravity + Kd*V[i] + F_wind;
}
```

## Adding collision detection ##
Adding collision detection is a bit more involved and there are solid methods in literature for this. For example, Dr. Eishcen has gracefully shared his chapter on cloth collision detection from the book Cloth Modelling & Animation onlne here (http://legacy.mae.ncsu.edu/directories/faculty/eischen/AKPeters-Chap8.pdf) which contains a pretty fast collision function this may help as a starting point. I will see if I could do a simple demo on this. But for the time being, this should get the interested individuals get started.