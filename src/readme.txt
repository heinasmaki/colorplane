Compile the code by issuing "make".

The code requires a fortran compiler. The makefile defaults to gfortran, which is recommended. 
The code has been tested with GCC versions 5.4 and 7.4. The choice of the fortran compiler is important as it affects the 
pseudorandom number generator (prng). The code uses the fortran built-in prng, which in GCC 5.4 is the Mersenne twister, 
and in GCC 7.4 xorshift1024*. Both generators are of high quality and give consistent results, but if you use some other 
fortran compiler, beware of the prng issue. 

All the results found in the "colorplane/materials" subfolder have been obtained with GCC 5.4.

The folder includes two auxiliary (GNU/Octave or Matlab) routines. They are used to compute the relative energies for lattice sites,
which can then be displayed as heat maps, similar to what is seen in the paper. They also produce the energy curves,
as seen in Fig. 5 of the paper. 

You can also use Octave or Matlab for displaying the lattices (and heat maps). For a lattice file "mylat.lat" issue
L = load ("mylat.lat");
imagesc (L)
