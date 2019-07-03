!
!   Colorplane: simulated annealing of a modified Potts model
!
!   Copyright (C) 2018 Sami Heinäsmäki
!
!   Colorplane is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Colorplane is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Colorplane.  If not, see <https://www.gnu.org/licenses/>.
!
program colorplane
  use colorplane_consts
  use colorplane_utils
  implicit none

  type(basedata) :: lattice
  integer(kind=ip), dimension(:,:), allocatable :: spins
  integer :: i, j, k, ct, m, cx, cy, id, is, js, nloops, loopind
  real(kind=dp) :: delsum, T, e_old, e_new, r, p, err
  integer(kind=ip) :: s_old, s_new, s_ref, s_out
  real :: start, finish  ! for timing
  integer, parameter :: latfil = 33, errfil = 99
  real(kind=dp), dimension(:), allocatable :: cumerr

  call collect_input (lattice)
  call check_if_restart (lattice)
  call finalize_input (lattice)
  
  ! Initialize spins: an outer padding of zeros (never updated)
  ! which is needed for the fixed-range energy calculation to work,
  ! and the inner "actual" Ny x Nx lattice of random integers 1...q
  ! either randomly or from a previous file
  allocate (spins(lattice%Ny + 2*lattice%d, lattice%Nx + 2*lattice%d))
  spins = 0
  call create_initial_spins ( lattice, spins )

  ! Finally, if we are doing a restart, change output name:
  if (lattice%restart) lattice%outfil = trim(lattice%outfil)//'-res'
    
  ! the average number of dist=1 equal spins at every step
  ! can be helpful in determining if the number of constant T
  ! update loops is high enough:
  nloops = updatenum (lattice%tmax, lattice%ctmax, lattice%offset)
  allocate( cumerr(nloops) )

  ! Start the actual annealing sequence, traverse the spin lattice in
  ! typewriter order
  call cpu_time ( start )
  k       = 0
  delsum  = sum (lattice%dval)
  loopind = 0
  
  do while (k < lattice%tmax)

     k = k + 1
     T = lattice%tvec(k)

     ! constant temperature MC loop
     ct = 0

     ! perform more constant T loops in the beginning, where
     ! the temperature is high, gradually shorten the below
     ! loop according to parameter offset
     do while (ct < (lattice%ctmax - lattice%offset * k) )

        loopind = loopind + 1
        ct      = ct + 1
        err     = 0.0_dp  ! average site error
        m       = 0       ! site counter
   
        ! lattice update loops
        do i = 1, lattice%Nx
           do j = 1, lattice%Ny

              ! We look at spin at the center of a subarray,
              ! whose upper-left corner is at (j,i).
              ! x and y coordinates of the center spin:
              cx    = i + lattice%d
              cy    = j + lattice%d
              s_old = spins(cy, cx)

              ! randomly generate s_new /= s_old:
              s_new = newspin (s_old, lattice%q)
              ! for q>~6 below is marginally faster:
              !call random_number (r)
              !ind   = 1 + floor( (lattice%q - 1) * r )
              !s_new = lattice%spinmat(s_old, ind)

              ! Compute site energies by looping over the nonzero
              ! values of the delta function in the sublattice
              ! around the site (cy, cx)          
              m     = m + 1
              e_old = 0.0_dp
              e_new = 0.0_dp
              do id = 1, lattice%dnum
                 is    = i + lattice%xind(id) - 1
                 js    = j + lattice%yind(id) - 1
                 s_ref = spins (js, is)
                 if (s_ref == s_old) e_old = e_old + lattice%dval(id)
                 if (s_ref == s_new) e_new = e_new + lattice%dval(id)
              end do

              ! Metropolis update of the spin
              s_out = s_old
              if (e_new < e_old) then
                 s_out = s_new
              else
                 p = exp(-(e_new - e_old) / T)
                 call random_number (r)
                 if (r < p) s_out = s_new
              end if

              spins(cy, cx) = s_out

              ! Cumulative average energy per lattice site
              err = err * (m-1)/real(m) + e_old/real(m)

           end do
        end do

        ! dist=1 equal spin values = normalized energies
        cumerr(loopind) = err / delsum

     end do

     ! print out ratio of parallel spins at dist = 1 at the end of
     ! constant T update cycle
     lattice%errvec(k) = cumerr(loopind)
     write(*,'(a,i4,a,F15.12)') 'Average number of neighbours per site' // &
          ' after step', k, ':', lattice%errvec(k)

  end do

  call cpu_time ( finish )
  write(*,'(/,a,F12.3,a)') 'Annealing time =', finish-start, ' seconds'

  call write_summary (lattice, k, finish-start)

  ! write out also the whole lattice, for creating a nice picture
  ! use e.g. Octave: B = load ("outfil.lat"); imagesc (B)
  open(latfil, file=trim(lattice%outfil)//'.lat', status='replace', & 
       action='write', form='formatted')
  
  do j = 1 + lattice%d, lattice%Ny + lattice%d
     write(latfil,*) (spins(j,i), i = 1+lattice%d, lattice%Nx + lattice%d)
  end do
  close (latfil)

  ! If requested, write out numbers of dist=1 equals in a big file
  if (lattice%writeneigh) call print_neighs (spins, lattice)

  ! cumerr is also written into a big file:
  open(errfil, file=trim(lattice%outfil)//'.dat', status='replace', & 
       action='write', form='formatted')
  do j = 1, nloops
     write(errfil, *) cumerr(j)
  end do
  close(errfil)
  
end program







