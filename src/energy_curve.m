%%
%%   Copyright (C) 2018 Sami Heinäsmäki
%%   E-mail: sami.heinasmaki@gmail.com
%%
%%   This is free software: you can redistribute it and/or modify
%%   it under the terms of the GNU General Public License as published by
%%   the Free Software Foundation, either version 3 of the License, or
%%   (at your option) any later version.
%%
%%   The code is distributed in the hope that it will be useful,
%%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%%   GNU General Public License for more details.
%%
%%   You should have received a copy of the GNU General Public License
%%   along with this code.  If not, see <https://www.gnu.org/licenses/>.
%%
function [y, M] = energy_curve (lattice, a)
  %%
  %% [y, M] = energy_curve (lattice, a)
  %% computes a curve y showing the average number of similar
  %% spins distance one away from the sites of the input lattice and
  %% for lattice constant a, as a function of distance from the center
  %% of the lattice. It is a cumulative average of the neighbour
  %% number matrix from inside out.
  %% Also outputs the neighbour matrix m if requested.

  %% size can also be non-square
  [nr, nc] = size (lattice);
  np       = floor (min([nr/2 nc/2]));
  y        = zeros (1, np);

  %% first use energy_map to create the matrix of averages
  M = energy_map (lattice, a);

  %% average M cumulatively from inside out:
  for p = 1:np

    a = nr-p+1;
    b = nc-p+1;
    [r, c] = size (M(p:a, p:b));
    y(np-p+1) = sum (sum (M(p:a, p:b))) / (r*c);

  end
  
endfunction

  

  

