function phi = lattice_delta(r)
  % Compact cosine kernel, matching the default Fortran interaction.
  phi = zeros(size(r));
  inside = abs(r) < 2;
  phi(inside) = 0.25 * (1 + cos(pi * r(inside) / 2));
end
