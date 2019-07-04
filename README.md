# Colorplane
Spin model simulator for plane colorings.

The principles and methods used by the code are described in the associated paper "Hamiltonian approach to plane colorings" in a file "plane_coloring.pdf" in the subfolder "materials".

The code is interactive, and you can look at the input file examples in the "materials/inputs" subfolder for reasonable parameters. 

The source code in the folder "src" requires a fortran compiler. The makefile defaults to gfortran, which is recommended. The code has been tested with GCC versions 5.4 and 7.4. The choice of the fortran compiler is important as it affects the pseudorandom number generator (prng). The code uses the fortran built-in prng, which in GCC 5.4 is the Mersenne twister, and in GCC 7.4 xorshift1024*. Both generators are of high quality and give consistent results, but if you use some other fortran compiler, beware of the prng issue. All the results found in the "materials" subfolder have been obtained with GCC 5.4.

The source code in the folder "src" is subject to the GPL 3.0 License. See the license file in the main folder for the full text.

The associated materials in the folder "materials" are licensed under a
Creative Commons Attribution-ShareAlike 4.0 International License.

These include especially the file "plane_coloring.pdf" and the lattices in the subfolder "lattices".

See the file "materials/License_CC" for the full text of the CC license.
