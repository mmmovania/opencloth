# Introduction #

This text will try to identify differences between the Jackobsen's approach and position based dynamics.


# Explanation #
Basically, for any method you use for cloth simulation, it has three things that co-exist in the following order:
  1. The calculationg of external forces
  1. The application of constraints (i.e. calculation of internal forces)
  1. Integration of positions and velocities
The simulation tries to solve the 2nd order ODE (Ma+Cv+Kx=F\_ext) using an integration of 2 first order differential eq. One for position and other for velocity. So if you use explicit Euler integration the equations are:
```
v(t+dt) = v(t) + dt*a(t)
x(t+dt) = x(t) + dt*v(t)
```
Verlet is a method for integration that approximates the simulation ODE by using a single equation for calculating the new position and implicitly represents the velocity using:
```
x(t+dt) = 2x(t) - x(t-dt) + a(t)dt^2
```
Now posiiton based dynamics(PBD) is a method of solving the 2nd order equation that works by first calculating the initial guesses of new positions (X\_pred) using explicit integration. Since explicit integration is not accurate enough, therefore, it applies constraints (distance/bending/collision constraints etc.) on these predicted positions to make sure that the predicted values are valid. Once the predicted values pass these constraints, they are integrated to the new position. Details about this method are given in the excellent SIGGRAPH 2008 coursenotes (Realtime Physics: http://www.matthiasmueller.info/realtimephysics/)

Now we can use verlet or any other integration method in the final step of PBD as well to get the new position.
So the difference is clear verlet is an integration scheme, position based dynamics is an algorithm that may use verlet integration in the final step to get the new positions.

I hope this clears up the differences between the two. <br>
Additional read: <br>
Jackobsen's approach(<a href='http://www.gamasutra.com/resource_guide/20030121/jacobson_01.shtml'>http://www.gamasutra.com/resource_guide/20030121/jacobson_01.shtml</a>)<br>
PBD (<a href='http://www.matthiasmueller.info/publications/posBasedDyn.pdf'>http://www.matthiasmueller.info/publications/posBasedDyn.pdf</a>)<br>