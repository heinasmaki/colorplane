Lattices are ascii files of integers. Most easily processed with e.g. GNU octave via 

%% Loading the lattice into a matrix L:
L = load ("filename.lat");

%% Showing a colored picture of the lattice:
imagesc (L)

Python alternatives (from the repository root):
python -m pip install -r scripts/requirements-plot.txt
python scripts/plot_lattice.py materials/lattices/q4-10-4-res-res.lat --spacing .0125 --output figures/lattice.png
python scripts/plot_neighbors.py materials/lattices/q4-10-4-res-res.lat --spacing .0125 --output figures/neighbors.png

See the main README for probability normalization, other spacings, and plotting
options. Use the actual spacing associated with the chosen lattice.
