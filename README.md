# Colorplane

Spin model simulator for plane colorings. Computes the ground state of a multicomponent spin system in a lattice, where the number of equal spin values unit distance apart is minimized. When interpreting spin values as colors, the output is an approximation to a [plane coloring](https://en.wikipedia.org/wiki/Hadwiger%E2%80%93Nelson_problem). 

The principles and methods used by the code are described in the associated paper [Hamiltonian approach to plane colorings](https://arxiv.org/abs/1908.03880).

The code is interactive, and you can look at the input file examples in the "materials/inputs" subfolder for reasonable parameters. 

The source code in the folder "src" requires a fortran compiler. The makefile defaults to gfortran, which is recommended. The code has been tested with GCC versions 5.4 and 7.4. The choice of the fortran compiler is important as it affects the pseudorandom number generator (prng). The code uses the fortran built-in prng, which in GCC 5.4 is the Mersenne twister, and in GCC 7.4 xorshift1024*. Both generators are of high quality and give consistent results, but if you use some other fortran compiler, beware of the prng issue. All the results found in the "materials" subfolder have been obtained with GCC 5.4.

You should think this software as a reference implementation: I would actually like you to write your own simulator which could use different choices of the simulated annealing algorithm. This would provide more information and could even give some important insight into the problem.

The source code in the folder "src" is subject to the GPL 3.0 License. See the license file in the main folder for the full text.

The associated materials in the folder "materials" are licensed under a
Creative Commons Attribution-ShareAlike 4.0 International License.

These include especially the paper and the lattices in the subfolder "lattices".

You can access [the full text](materials/License_CC) of the CC license.
