!
!   Copyright (C) 2018 Sami Heinäsmäki
!
!   This file is part of Colorplane.
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
module colorplane_consts
  implicit none
 
  ! For the spin lattice we use the shortest possible integers.
  ! One-byte integers are good for spin values 1 - 255
  ! (zero value is reserved for other use). 
  ! If this is not enough, modify here the accuracy ip
  integer, parameter :: ip = selected_int_kind (1)  

  ! sp is the standard real, typically 32 bit
  integer, parameter :: sp  = kind(1.0)
  ! and dp is the "long real", usually 64 bit
  integer, parameter :: dp  = selected_real_kind(2*precision(1.0_sp))

  ! Some constants for the mc computation:

  ! Critical value of relative difference of old and new lattice
  ! energies to stop annealing:
  real(kind=dp), parameter :: err_lim = 1.0e-12

  ! Minimum value of delta function at the lattice site
  real(kind=dp), parameter :: del_min = 1.0e-14

  ! Parameter offset reduces the number of constant T loops as the
  ! temperature is lowered, according to rule
  ! N_ctloops = Nctmax - offset * k
  ! where k is the annealing step.
  ! This can be either zero or some small integer, experience
  ! has shown that convergence takes longer time at high T.
  ! This can optionally be set at initialization.
  integer, parameter, public :: offset_default = 0

  ! Datatype basedata holds the basic information about the problem. 
  ! Defaults are overwritten once input data is processed.  
  type :: basedata
     ! number of spins = colors:
     integer(kind=ip) :: q
     
     ! lattice dimensions (default is to use Nx = Ny):
     integer :: Nx, Ny
     
     ! lattice spacing:
     real(kind=dp) :: a
     
     ! max number of annealing steps
     integer :: tmax
     
     ! reduce temperature at each annealing step by 0 < tscale < 1:
     real(kind=dp) :: tscale
     
     ! number of update steps at constant temperature:
     integer :: ctmax
     
     ! root name for output files; there will be outfil.sum for basic data,
     ! outfil.lat for the spin lattice itself, and optionally outfil.err for
     ! errors site by site
     character(len=256) :: outfil
     
     ! number giving the maximum number of lattice points where the
     ! delta function interaction extends from any site
     integer :: d

     ! vector to hold the values of the lattice delta function,
     ! giving the interaction strength with the central spin
     real(kind=dp), dimension(:), allocatable :: dval

     ! vectors to hold the x and y indices of the nonzero values
     ! of the delta function, in the (2d+1)x(2d+1) subarray
     integer, dimension(:), allocatable :: xind, yind

     ! length of the previous vectors
     integer :: dnum

     ! if restart is false, do a fresh computation based on input,
     ! otherwise a restart using previous output as a starting point
     logical :: restart = .false.

     ! vector to hold the temperature values,
     ! one value for each annealing step
     real(kind=dp), dimension(:), allocatable :: tvec

     ! minimum temperature in annealing
     real(kind=dp) :: t_min
     
     ! vector to hold the error = average number of parallel spins
     ! at distance 1 away, one value for each annealing step
     real(kind=dp), dimension(:), allocatable :: errvec

     ! q x (q-1) array to hold spin values where a new spin value
     ! can be drawn during lattice update steps
     integer(kind=ip), dimension(:,:), allocatable :: spinmat

     ! whether to write numbers of neighbours for all sites out
     logical :: writeneigh = .false.

     ! offset value, can optionally be set to non default value
     integer :: offset = offset_default
     
  end type basedata

end module colorplane_consts
