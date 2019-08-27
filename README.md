# Colorplane

Spin model simulator for plane colorings. Computes the ground state of a multicomponent spin system in a lattice, where the number of equal spin values unit distance apart is minimized. When interpreting spin values as colors, the output is an approximation to a [plane coloring](https://en.wikipedia.org/wiki/Hadwiger%E2%80%93Nelson_problem). 

The principles and methods used by the code are described in the associated paper [Lattice approach to plane colorings](https://arxiv.org/abs/1908.03880).

The code is interactive, and you can look at the input file examples in the "materials/inputs" subfolder for reasonable parameters. 

You should regard this software as a reference implementation: I would actually like you to write your own simulator which could use different choices of the simulated annealing algorithm. This would provide more information and could even give some important insight into the problem.

The source code in the folder "src" is subject to the GPL 3.0 License. See the license file in the main folder for the full text. The associated files in the folder "materials" and its subfolders are licensed under a Creative Commons Attribution-ShareAlike 4.0 International License. The material includes the input files used in the paper, as well as the lattices shown there. You can access [the full text](materials/License_CC) of the CC license.
