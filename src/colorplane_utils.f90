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
module colorplane_utils

  use colorplane_consts
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none

  public :: collect_input
  ! interactively collect input data from user

  public :: check_if_restart
  ! looks for previous output, if found, assumes a restart

  public :: finalize_input
  ! precomputes data for the mc iteration based on input

  public :: create_initial_spins
  ! reads spins from previous output, or creates a new random set

  private :: lattice_delta
  ! lattice delta function, approximation based on cosine

  private :: set_annealing_scheme
  ! determines the temperature values in annealing

  public :: newspin
  ! returns a new spin different from old spin, based on random_number

  public :: updatenum
  ! computes the total number of lattice updates

  public :: print_neighs
  ! computes the the errors for all lattice sites and prints them

  public :: write_summary
  ! Prints out a summary: temperature and probability that spins
  ! have equal values at distance 1 away from each other


contains

  subroutine collect_input (lattice_data)
    !
    type(basedata), intent(out) :: lattice_data
    character(len=1) y_or_n
    logical :: not_complete
    integer :: ierror, ncols

    write(*,'(/,a,/)') '             Welcome to COLORPLANE'
    write(*,'(a)') 'Simulated annealing of a fixed-range' // &
         ' interaction Potts model'
    write(*,'(a)') 'Copyright (C) 2018 Sami Heinäsmäki'
    write(*,'(a,/)') 'Colorplane is free software, see ' // &
    '<https://www.gnu.org/licenses/>'

    not_complete = .true.
    do while (not_complete)
       write(*,'(/,a)') 'Give the number of colors, from 2 to 255'
       read(*,'(i3)',iostat=ierror) lattice_data%q
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read integer value, please redo...'
          cycle
       end if
       if (lattice_data%q > 1 .and. lattice_data%q <= 255) not_complete = .false.
    end do

    not_complete = .true.
    do while (not_complete)
       write(*,'(/,a)') 'Enter number of rows in a square lattice'
       read(*,'(i8)',iostat=ierror) lattice_data%Ny
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read integer value, please redo...'
          cycle
       end if
       if (lattice_data%Ny > 1) not_complete = .false.
    end do
    ! default is square lattice:
    lattice_data%Nx = lattice_data%Ny

    not_complete = .true.
    do while (not_complete)
       write(*,'(/,a)') 'Enter lattice spacing'
       read(*,*,iostat=ierror) lattice_data%a
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read real value, please redo...'
          cycle
       end if
       if (lattice_data%a > 0.0) not_complete = .false.
    end do

    not_complete = .true.
    do while (not_complete)
       write(*,'(/,a)') 'Enter max number of annealing steps'
       read(*,'(i6)',iostat=ierror) lattice_data%tmax
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read integer value, please redo...'
          cycle
       end if
       if (lattice_data%tmax > 0) not_complete = .false.
    end do

    not_complete = .true.
    do while (not_complete)
       write(*,'(/,a)') 'Enter temperature scaling between the'// &
            ' annealing steps, in interval (0,1)'
       read(*,*,iostat=ierror) lattice_data%tscale
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read real value, please redo...'
          cycle
       end if
       if (lattice_data%tscale > 0.0 .and. lattice_data%tscale < 1.0) then
          not_complete = .false.
       end if
    end do

    not_complete = .true.
    do while (not_complete)
       write(*,'(/,a)') 'Enter number of update steps at constant temperature'
       read(*,'(i6)',iostat=ierror) lattice_data%ctmax
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read integer value, please redo...'
          cycle
       end if
       if (lattice_data%ctmax > 0) not_complete = .false.
    end do

    write(*,'(/,a)') 'Enter a root name for output data;'
    write(*,'(a)')   'file with extension .sum contains basic output data,'
    write(*,'(a)')  'file with extension .lat contains the final lattice,'
    write(*,'(a)')  'and the optional file with extension .ngb contains'// &
         ' neighbours per site.'
    read(*,'(a)') lattice_data%outfil

    write(*,'(/,a)') 'Print out number of neighbours for all sites? (y/N)'
    write(*,'(a)') 'Warning: output may be large!'
    read(*,'(a1)') y_or_n
    if (y_or_n == 'y' .or. y_or_n == 'Y') lattice_data%writeneigh = .true.

    write(*,'(/,a)') 'Modify defaults? (y/n)'
    read(*,'(a1)') y_or_n
    if (y_or_n == 'y') then
       write(*,'(/,a)') 'For a non-square lattice,'// &
            ' enter the number of columns in lattice, or zero'// &
            ' to stick with default = number of rows .'
       read(*,*) ncols
       if (ncols > 1) lattice_data%Nx = ncols
       ! offset inquiry
       write(*,'(/,a,i2,a)') 'The default is to reduce the number'// &
            ' of constant temperature updates by ', offset_default, &
            ' at each time temperature is reduced, set new value'// &
            ' or empty to stick with default:'
       read(*,'(i4)',iostat=ierror) lattice_data%offset
       if (ierror < 0) error stop 'Unexpected end of input'
       if (ierror /= 0) then
          write(*,*) 'Could not read integer value, using default value.'
          lattice_data%offset = offset_default
       end if
    end if

    write(*,*)

  end subroutine collect_input


  subroutine check_if_restart (data)
    !
    ! Checks if files root.lat and root.sum are present,
    ! assumes a restart computation in that case
    !
    type(basedata), intent(inout) :: data
    logical :: lat_exists, sum_exists

    inquire (file=trim(data%outfil)//'.lat', exist=lat_exists)
    inquire (file=trim(data%outfil)//'.sum', exist=sum_exists)

    if (lat_exists .and. sum_exists) then

       data%restart = .true.

       write(*,'(/,a)') 'Previous output found, assuming restart...'

    end if

  end subroutine check_if_restart


  subroutine finalize_input (data)
    !
    ! Precomputes the delta function data based on input,
    ! also fixes temperature scheme
    !
    type(basedata), intent(inout) :: data
    integer :: m, dd, i, j, dc
    real(kind=dp) :: dist
    integer, dimension(:), allocatable :: xtmp, ytmp
    real(kind=dp), dimension(:), allocatable :: dtmp

    ! Extent of delta interaction from a given lattice site
    ! The lattice delta function interaction ranges so that
    ! -2 < |xi - xj| -1 < 2, and when |xi - xj| = ad,
    ! we get d = ceiling (2 + 1/a)
    data%d = ceiling (2.0_dp + 1.0_dp / data%a)

    ! m^2 is the upper limit for amount of data values
    m = 2*data%d + 1
    allocate ( dtmp(m**2), xtmp(m**2), ytmp(m**2) )

    ! Precompute the lattice delta function.
    ! Note that the delta function is purely geometrical (independent
    ! of the lattice spins) and the nonzero values and their locations
    ! can be computed once and for all for a sublattice, where the
    ! site spin sits in the middle.
    dd = data%d + 1
    dc = 0
    do i = 1, m
       do j = 1, m
          dist = sqrt (real(dd-i,dp)**2 + real(dd-j,dp)**2)
          dist = lattice_delta(dist - 1.0_dp/data%a, data%kernel)
          dist = dist / data%a
          if (dist > del_min) then
             dc       = dc + 1
             xtmp(dc) = i
             ytmp(dc) = j
             dtmp(dc) = dist
          end if
       end do
    end do

    ! delta function values and positions go in struct data:
    data%dnum = dc
    allocate ( data%dval(dc), data%xind(dc), data%yind(dc) )
    data%dval = dtmp(1:dc)
    data%xind = xtmp(1:dc)
    data%yind = ytmp(1:dc)

    ! errvec holds the numbers of parallel-spin neighbours
    ! at unit distance for each temperature step
    allocate( data%errvec(data%tmax) )
    data%errvec = -1.0_dp

    ! The choice of the temperature values is given in its
    ! own subroutine; modify it if you feel to
    if (dc == 0) error stop 'Empty interaction stencil'
    call set_annealing_scheme (data)

  end subroutine finalize_input


  subroutine create_initial_spins (data, SM)
    !
    ! either creates a random integer matrix (new run) or
    ! reads the matrix of spin values based on previous output
    !
    type(basedata), intent(in) :: data
    integer(kind=ip), dimension(:,:), intent(inout) :: SM
    integer, parameter :: oldlat = 51
    integer :: i, j
    real(kind=dp) :: r

    if (data%restart) then

       open(oldlat, file=trim(data%outfil)//'.lat', status='old', &
            action='read', form='formatted')

       do j = 1 + data%d, data%Ny + data%d
          read(oldlat,*) (SM(j, data%d + i), i = 1, data%Nx)
       end do

       close (oldlat)

       write (*,'(/,a)') 'old spins read from file...'

    else

       do i = 1 + data%d, data%Nx + data%d
          do j = 1 + data%d, data%Ny + data%d
             call random_number (r)
             SM(j, i) = int(1 + floor(data%q * r),ip)
          end do
       end do

       write (*,'(/,a,/)') 'new random spins created...'

  end if

  end subroutine create_initial_spins


  function lattice_delta(r, kernel) result(phi)
    !
    ! Returns approximate for the lattice delta function
    !
    real(kind=dp), intent(in) :: r
    character(len=*), intent(in) :: kernel
    real(kind=dp) :: pi, phi

    pi  = 3.1415926535897932_dp
    phi = 0.0_dp

    if (kernel == 'hat') then
       phi = max(0.0_dp, 1.0_dp - abs(r))
    else if (abs(r) < 2.0_dp) then
       phi = .25_dp * (1 + cos (pi * r / 2.0_dp))
    end if

  end function lattice_delta


  subroutine set_annealing_scheme (data)
    !
    ! Set the initial temperature and later temperature values
    ! based on user input.
    !
    type(basedata), intent(inout) :: data
    real(kind=dp) :: prefac
    integer :: i

    allocate( data%tvec(data%tmax) )

    ! In the thermally dominated case we expect that each site has
    ! on average 1/q parallel-spin neighbours.
    ! The energy in this situation can be computed by taking the
    ! sum over delta function values and dividing by q.
    ! Experimenting has shown that this may be higher than
    ! necessary: the simulation does not converge in the beginning.
    ! Therefore, this initial value is reduced by multiplying
    ! with a prefactor < 1
    prefac = 0.1_dp

    ! A simple scaling scheme is employed here, there are
    ! alternatives in the literature but this is easy for
    ! the user to affect, as the scaling factor 0 < tscale < 1
    ! is queried at input
    if (data%restart) then

       ! Full-precision checkpoint replaces parsing rounded human output.
       data%tvec(1) = data%tscale * data%restart_temperature

       write (*,'(/,a,/)') 'initial temperature read from file...'

    else

       data%tvec(1) = prefac * sum(data%dval) / data%q

    end if

    ! Finally we have the initial temperature,
    ! scale down from that
    do i = 2, data%tmax
       data%tvec(i) = data%tscale * data%tvec(i-1)
    end do
    ! Use global step indices so split runs have the same schedule.
    if (data%reheat_every > 0) then
       do i = 1, data%tmax
          if (mod(data%completed_steps+i-1, data%reheat_every) == 0 .and. &
               data%completed_steps+i > 1) data%tvec(i:) = data%tvec(i:) * data%reheat_factor
       end do
    end if
    if (any(data%tvec <= tiny(1.0_dp)) .or. any(.not. ieee_is_finite(data%tvec))) &
         error stop 'Temperature out of range: shorten cooling schedule'

  end subroutine set_annealing_scheme


  function newspin (si, q) result(sf)
    !
    ! Returns sf /= si,  1 <= sf <= q
    !
    integer(kind=ip), intent(in) :: si, q
    integer(kind=ip) :: sf
    real(kind=dp) :: r
    call random_number(r)
    sf = int(1 + floor(real(q-1,dp)*r), ip)
    if (sf >= si) sf = sf + 1_ip

  end function newspin


  function updatenum (Na, Nct, ofs) result (nloops)
    !
    ! Returns the total number of lattice updates
    !
    integer, intent(in) :: Na, Nct, ofs
    integer :: nloops
    integer(kind=selected_int_kind(18)) :: steps,total
    steps=int(Na,kind(steps))
    if (Na<1 .or. Nct<1 .or. ofs<0) error stop 'Invalid sweep schedule'
    if (int(Nct,kind(steps))-steps*ofs<1) error stop 'Empty temperature stage'
    total=steps*Nct - (steps*(steps+1)/2)*ofs
    if (total>huge(nloops)) error stop 'Sweep count overflow'
    nloops=int(total)

  end function updatenum


  subroutine print_neighs (spins, data)
    !
    ! Prints the numbers of parallel-spin neighbours separated
    ! by unit distance for all lattice sites into a file
    !
    integer(kind=ip), dimension(:,:), intent(in) :: spins
    type(basedata), intent(in) :: data

    integer(kind=ip) :: sc, ss
    integer :: i, j, cx, cy, is, js, k, id
    real(kind=dp) :: norm, err
    real(kind=dp), dimension(data%Nx) :: nvec
    integer, parameter :: chan = 55

    open(chan, file=trim(data%outfil)//'.ngb', status='replace', &
       action='write', form='formatted')

    ! Normalization factor is the sum of all
    ! delta function values
    norm = sum (data%dval)

    do j = 1, data%Ny
       do i = 1, data%Nx

          ! x and y coordinates of the center spin of
          ! the sublattice
          cx = i + data%d
          cy = j + data%d

          ! spin at the center of a ns x ns sublattice,
          ! whose upper-left corner is at (j,i)
          sc = spins(cy, cx)

          ! Compute site energy: loop over nonzero delta
          ! function values in the sublattice around sc
          err = 0.0_dp
          do id = 1, data%dnum
             is = i + data%xind(id) - 1
             js = j + data%yind(id) - 1
             ss = spins (js, is)
             if (ss == sc) err = err + data%dval(id)
          end do

          nvec(i) = err / norm

       end do

       ! write one line
       write (chan, *) (nvec(k), k = 1, data%Nx)

    end do

    close (chan)

  end subroutine print_neighs


  subroutine write_summary (data, k, elapsed_time)
    !
    ! write out basic data to .sum file
    !
    type(basedata), intent(in) :: data
    integer, intent(in) :: k
    real, intent(in) :: elapsed_time

    integer, parameter :: output = 22
    character(len=8) :: date
    character(len=10) :: time
    integer :: i
    real(kind=dp) :: xsize, ysize

    ! compute also the size of the lattice in unit distances,
    ! because this is essential for the coloring problem
    xsize = data%a * real(data%Nx,dp)
    ysize = data%a * real(data%Ny,dp)

    call date_and_time(date,time)

    ! Open the output file and print the header:

    open(output, file=trim(data%outfil)//'.sum', status='replace', &
         action='write', form='formatted')

    write(output,'(54a1)') ('=', i=1,54)
    write(output,*) 'Summary of the COLORPLANE output:'
    write(output,'(54a1,/)') ('=', i=1,54)

    write(output,*) 'COLORPLANE run on ', date(1:4),'-',date(5:6),'-', &
         date(7:8), ' at ', time(1:2),':',time(3:4)

    write(output,'(/,a,/)') 'The following input data was used:'

    write(output,'(a,i0)')      'Number of spins (=colors):      ', data%q

    write(output,'(a,i4,a,i4,a)') 'Lattice size:                   ', &
         data%Ny, ' x', data%Nx, ' lattice sites'

    write(output,'(a,f6.2,a,f6.2,a)') '                              = ', &
         ysize, ' x', xsize, ' unit distances'

    write(output,'(a,f8.5)')    'Lattice spacing:                ', data%a

    write(output,'(a,i3)')      'Number of annealing loops:      ', k

    write(output,'(a,i4)')      'Number of updates at constant T:', data%ctmax

    write(output,'(a,f6.3)')   'Temperature scaling:            ', data%tscale

    write(output,'(a,f12.2,/)') 'Elapsed time (seconds):     ', elapsed_time

    ! Write out the results by annealing step:
    write(output,'(54A1)') ('-', i=1,54)
    write(output,'(a)') 'Average number of parallel-spin neighbours,'
    write(output,'(a)') 'separated by unit distance:'
    write(output,'(54A1,/)') ('-', i=1,54)

    write(output,'(a7,a16,a18,/)') 'step', 'temperature', 'parallel spins'

    do i=1,k
       write(output,'(i8,2es25.17)') i, data%tvec(i), data%errvec(i)
    end do

    close(output)

  end subroutine write_summary

end module colorplane_utils
