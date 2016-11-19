# apc_hw4

This assignment solves the heat diffusion equation partialdvT/partialdvt = kappa*delta^2*T 
with kappa = constant on a 2D domain of size 0 leq x leq pi and 0 leq y leq pi. The boundary
conditions are T(x,0) = cos^2(x), T(x,pi)=sin^2(x), T(0,y)=T(pi,y) (periodic in x).

The equation is solved using centered finite differences in space and the forward Euler 
method in time.

There are 3 different implementations provided:

heat_serial.cpp is run ./heat_serial <nx> for grid size nx^2. This is a serial implementation 
that solves the problem starting from i.c. T=0 everywhere. Contours or an image of the final
temperature and the final, volume averaged temperature are provided using grid sizes of 128^2,
256^2, and 512^2, where the solution is run until time t=0.5*pi^2/kappa. 

heat_omp is run ./heat_omp <nx> <nthreads> for grid size nx^2. This is the serial implementation updated using openMP. Same details as above.

heat_mpi is still being finished.