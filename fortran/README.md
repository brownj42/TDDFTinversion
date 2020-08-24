# Fortran usage

The main program is in tddftinversion.f90. 

To compile, you may have to change the LFLAGS line to point to the correct locations for the LAPACK, ARPACK and BLAS libraries then type make.

To run ./tddftinversion.

This is setup such that you can change the [parameters.in](parameters.in) file to change the system parameters.

For the two particle with driving potential example in the paper
<pre>
1              number of dimensions each particle live in  
71             number of one-dimensional grid points  
2              number of particles  
-14            min x-value of grid    
14             max x-value of grid    
0.005          time step  
0              0=Pseudo Inverse, 1=MINRES-QLP  
</pre>

For the 1 particle 3D example
<pre>
3              number of dimensions each particle live in  
25             number of one-dimensional grid points  
1              number of particles  
-6.            min x-value of grid  
6.             max x-value of grid  
0.00314159     time step  
1              0=Pseudo Inverse, 1=MINRES-QLP  
</pre>

For the 1 particle 1D splitting wavepacket example
<pre>
1              number of dimensions each particle live in  
111            number of one-dimensional grid points  
1              number of particles  
-11            min x-value of grid  
11             max x-value of grid  
0.00314159     time step  
0              0=Pseudo Inverse, 1=MINRES-QLP  
</pre>
