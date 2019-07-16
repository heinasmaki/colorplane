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
function M = energy_map (lattice, a)
  %%
  %% computes the matrix containing the average number of similar
  %% spins distance one away from the sites of the input lattice and
  %% for lattice constant a

  %% size can also be non-square
  [nr, nc] = size (lattice);
  M        = zeros(nr, nc);

  %% the computation needs the delta function:
  %% The lattice delta function interaction ranges so that -2 < |xi -
  %% xj| -1 < 2, and when |xi - xj| = ad, we get d = ceil(2 + 1/a)
  d = ceil (2 + 1 / a);
    
  %% Precompute the lattice delta function for the array of size zpad
  %% to be used in the interaction energy calculations.
  %% Note that the delta function is purely geometrical (independent
  %% of the lattice spins) and its values can be computed once and for
  %% all for a sublattice, where the site spin sits in the middle.
  del_cut = 1e-12;   % discard smaller delta values
  [xx,yy] = meshgrid(-d:d);
  distmat = sqrt (xx.^2 + yy.^2);
  deltmat = (1 / a) * lattice_delta (distmat - 1/a);
  id      = find (deltmat > del_cut);
  dv      = deltmat(id);
  dsum    = sum(dv);

  %% Embed spins into a larger lattice of zeros
  spins = zeros (nr + 2*d, nc + 2*d);
  spins(1+d:nr+d, 1+d:nc+d) = lattice;
  
  err = 0;  m = 0;

  for jx = d + 1:nc + d
    for jy = d + 1:nr + d

      m    = m + 1;
      sval = spins(jy, jx);

      %% Take a sublattice around current site, containing all spins
      %% which interact with the site
      spad = spins(jy-d:jy+d,jx-d:jx+d);
      %% delta function interaction between spins distance 1 apart
      sdel = spad(id);     % sites where delta function nonzero
      is   = find (sdel ~= sval);  % find spins different from s0
      sdel(is) = 0;                % and make them zero
      M(jy-d,jx-d) = (sdel' * dv) / (sval * dsum);   % sum over lattice spins
      
    end
  end
  
endfunction

  

  

