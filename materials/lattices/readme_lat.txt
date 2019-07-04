Lattices are ascii files of integers. Most easily processed with e.g. GNU octave via 

%% Loading the lattice into a matrix L:
L = load ("filename.lat");

%% Showing a colored picture of the lattice:
imagesc (L)
