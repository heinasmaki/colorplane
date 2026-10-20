! Reproducibility, checkpoint validation, and measurements of fixed configurations.
! Copyright (C) 2026 Sami Heinasmaki. GPL-3.0-or-later, as the rest of Colorplane.
module colorplane_reliability
  use colorplane_consts
  use colorplane_build_info
  use, intrinsic :: iso_fortran_env, only: int64, compiler_version, compiler_options
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: read_options, validate_input, initialize_rng, read_checkpoint, write_checkpoint
  public :: load_lattice, save_lattice, measure, shuffle_sites, usage, write_metadata
  public :: lattice_hash, check_output_root, save_rng

contains

  subroutine usage()
    print '(a)', 'Usage: colorplane [options] < input.inp'
    print '(a)', '--seed N                 Nonnegative integer seed (default 1)'
    print '(a)', '--initial FILE           Start from a supplied lattice with a fresh RNG'
    print '(a)', '--order typewriter|random  Visit each site once per sweep'
    print '(a)', '--kernel cosine|hat      Interaction kernel (default cosine)'
    print '(a)', '--reheat-every N          Reheat before every N+1, 2N+1, ... stage (0 disables)'
    print '(a)', '--reheat-factor F         Reheating multiplier, greater than 1'
    print '(a)', 'An existing output root requires a matching .state checkpoint for restart.'
  end subroutine

  subroutine read_options(data)
    type(basedata), intent(inout) :: data
    character(len=1024) :: arg, val
    integer :: k, ios
    k = 1
    do while (k <= command_argument_count())
       call get_command_argument(k, arg)
       if (k == command_argument_count()) error stop 'Missing option value'
       call get_command_argument(k+1, val)
       ios = 0
       select case(trim(arg))
       case('--seed')
          read(val,*,iostat=ios) data%seed
          data%seed_given = .true.
       case('--initial')
          data%initial_file = val
       case('--order')
          if (val /= 'random' .and. val /= 'typewriter') error stop 'Unknown sweep order'
          data%random_order = val == 'random'
       case('--kernel')
          if (val /= 'cosine' .and. val /= 'hat') error stop 'Unknown kernel'
          data%kernel = val(1:len(data%kernel))
       case('--reheat-every')
          read(val,*,iostat=ios) data%reheat_every
       case('--reheat-factor')
          read(val,*,iostat=ios) data%reheat_factor
       case default
          error stop 'Unknown command-line option; see --help'
       end select
       if (ios /= 0) error stop 'Invalid numeric option'
       k = k + 2
    end do
  end subroutine

  subroutine validate_input(data)
    type(basedata), intent(in) :: data
    integer(int64) :: last, loops
    if (data%q < 2 .or. data%q > 255) error stop 'Colors must be in 2..255'
    if (min(data%Nx,data%Ny) < 2) error stop 'Lattice dimensions must exceed 1'
    ! Coarse spacings >= 1/2 allow self-interactions, incompatible with spin-delta updates.
    if (.not. ieee_is_finite(data%a)) error stop 'Nonfinite lattice spacing'
    if (data%a <= 0 .or. data%a >= 0.5_dp) error stop 'Spacing must lie in (0, 0.5)'
    if (.not. ieee_is_finite(data%tscale)) error stop 'Nonfinite cooling factor'
    if (data%tscale <= 0 .or. data%tscale >= 1) error stop 'Cooling factor must lie in (0,1)'
    if (data%tmax < 1 .or. data%ctmax < 1 .or. data%offset < 0) error stop 'Invalid sweep schedule'
    last = int(data%ctmax,int64) - int(data%offset,int64)*(int(data%tmax,int64)+data%completed_steps)
    if (last < 1) error stop 'Each temperature stage must have at least one sweep'
    loops = int(data%tmax,int64)*data%ctmax
    if (loops > huge(1)) error stop 'Too many sweeps for indexing'
    if (int(data%Nx,int64)*data%Ny > huge(1)) error stop 'Too many sites for indexing'
    if (1.0_dp/data%a > 10000.0_dp) error stop 'Interaction stencil too large'
    if (data%seed < 0 .or. data%seed > 2147483646_int64) error stop 'Seed must lie in 0..2147483646'
    if (data%reheat_every < 0) error stop 'Negative reheat interval'
    if (.not. ieee_is_finite(data%reheat_factor)) error stop 'Nonfinite reheating factor'
    if (data%reheat_factor <= 1) error stop 'Reheating factor must exceed 1'
  end subroutine

  subroutine initialize_rng(data)
    type(basedata), intent(in) :: data
    integer, allocatable :: state(:)
    integer :: n, i
    integer(int64) :: x
    call random_seed(size=n)
    allocate(state(n))
    ! Modular expansion avoids overflow and supports every advertised seed.
    x = data%seed
    do i = 1,n
       x = modulo(48271_int64*x+1_int64,2147483647_int64)
       state(i) = int(x)
    end do
    call random_seed(put=state)
  end subroutine

  subroutine save_rng(data)
    type(basedata), intent(in) :: data
    integer :: u,n
    integer, allocatable :: state(:)
    call random_seed(size=n)
    allocate(state(n))
    call random_seed(get=state)
    open(newunit=u,file=trim(data%outfil)//'.rng-initial',status='replace')
    write(u,*) n
    write(u,*) state
    close(u)
  end subroutine

  subroutine state_ok(ios)
    integer, intent(in) :: ios
    if (ios /= 0) error stop 'Invalid or truncated checkpoint'
  end subroutine

  subroutine read_checkpoint(data)
    type(basedata), intent(inout) :: data
    integer :: u, ios, q, nx, ny, ctmax, offset, interval, n, expected
    logical :: exists, order
    real(dp) :: a, scale, factor
    integer, allocatable :: state(:)
    character(len=4096) :: line
    character(len=16) :: kernel
    inquire(file=trim(data%outfil)//'.state',exist=exists)
    if (.not. exists) error stop 'Restart needs .state; use --initial FILE for legacy lattices'
    if (data%seed_given .or. len_trim(data%initial_file)>0) error stop 'Seed/initial conflicts with restart'
    open(newunit=u,file=trim(data%outfil)//'.state',status='old')
    read(u,'(a)',iostat=ios) line
    call state_ok(ios)
    if (line /= 'COLORPLANE_STATE_1') error stop 'Unknown checkpoint format'
    read(u,'(a)',iostat=ios) line
    call state_ok(ios)
    if (line /= compiler_version()) error stop 'Checkpoint compiler mismatch'
    read(u,'(a)',iostat=ios) line
    call state_ok(ios)
    if (line /= compiler_options()) error stop 'Checkpoint compiler options mismatch'
    read(u,'(a)',iostat=ios) line
    call state_ok(ios)
    if (line /= build_identity) error stop 'Checkpoint source/build mismatch'
    read(u,*,iostat=ios) q,nx,ny,a
    call state_ok(ios)
    ! Exact comparisons are intentional: checkpoint text round-trips binary64.
    if (q/=data%q .or. nx/=data%Nx .or. ny/=data%Ny .or. a/=data%a) &
         error stop 'Checkpoint lattice parameters mismatch'
    read(u,*,iostat=ios) kernel,order,scale,ctmax,offset,interval,factor
    call state_ok(ios)
    if (kernel/=data%kernel .or. (order.neqv.data%random_order) .or. scale/=data%tscale .or. &
         ctmax/=data%ctmax .or. offset/=data%offset .or. interval/=data%reheat_every .or. &
         factor/=data%reheat_factor) error stop 'Checkpoint annealing parameters mismatch'
    read(u,*,iostat=ios) data%seed,data%completed_steps,data%restart_temperature
    call state_ok(ios)
    if (data%completed_steps < 1 .or. .not. ieee_is_finite(data%restart_temperature)) &
         error stop 'Invalid checkpoint temperature/step'
    if (data%restart_temperature <= 0) error stop 'Invalid checkpoint temperature'
    read(u,*,iostat=ios) n
    call state_ok(ios)
    call random_seed(size=expected)
    if (n/=expected) error stop 'Checkpoint RNG size mismatch'
    allocate(state(n))
    read(u,*,iostat=ios) state
    call state_ok(ios)
    read(u,*,iostat=ios) data%checkpoint_hash
    call state_ok(ios)
    close(u)
    call random_seed(put=state)
  end subroutine

  subroutine write_checkpoint(data,spins)
    type(basedata), intent(in) :: data
    integer(ip), intent(in) :: spins(:,:)
    integer :: u,n
    integer, allocatable :: state(:)
    call random_seed(size=n)
    allocate(state(n))
    call random_seed(get=state)
    open(newunit=u,file=trim(data%outfil)//'.state',status='replace')
    write(u,'(a)') 'COLORPLANE_STATE_1'
    write(u,'(a)') compiler_version()
    write(u,'(a)') compiler_options()
    write(u,'(a)') build_identity
    write(u,*) data%q,data%Nx,data%Ny,data%a
    write(u,*) trim(data%kernel),' ',data%random_order,data%tscale,data%ctmax,data%offset, &
         data%reheat_every,data%reheat_factor
    write(u,*) data%seed,data%completed_steps+data%tmax,data%tvec(data%tmax)
    write(u,*) n
    write(u,*) state
    write(u,*) lattice_hash(spins,data)
    close(u)
  end subroutine

  function lattice_hash(spins,data) result(hash)
    integer(ip), intent(in) :: spins(:,:)
    type(basedata), intent(in) :: data
    integer(int64) :: hash
    integer :: i,j
    hash=1
    do j=1,data%Ny
       do i=1,data%Nx
          hash=modulo(hash*65599_int64+spins(j+data%d,i+data%d),2147483647_int64)
       end do
    end do
  end function

  subroutine check_output_root(data)
    type(basedata), intent(in) :: data
    character(len=12), parameter :: extensions(*) = [character(len=12) :: &
         '.lat','.sum','.state','.meta','.metrics.csv','.dat','.ngb','.initial.lat','.rng-initial']
    logical :: exists
    integer :: i
    do i=1,size(extensions)
       inquire(file=trim(data%outfil)//trim(extensions(i)),exist=exists)
       if (exists) error stop 'Output destination already exists; choose a new root'
    end do
  end subroutine

  subroutine load_lattice(path,spins,data)
    character(len=*), intent(in) :: path
    type(basedata), intent(in) :: data
    integer(ip), intent(inout) :: spins(:,:)
    integer :: u,j,ios
    integer :: row(data%Nx+1)
    character(len=:), allocatable :: line
    ! Each row must have exactly Nx integers, without wrapping into the next row.
    allocate(character(len=8*(data%Nx+1)) :: line)
    open(newunit=u,file=path,status='old',iostat=ios)
    if (ios /= 0) error stop 'Cannot open starting lattice'
    do j=1,data%Ny
       read(u,'(a)',iostat=ios) line
       if (ios /= 0) error stop 'Starting lattice has too few rows'
       read(line,*,iostat=ios) row(1:data%Nx)
       if (ios /= 0) error stop 'Invalid lattice row'
       if (any(row(1:data%Nx)<1) .or. any(row(1:data%Nx)>data%q)) error stop 'Invalid spin value'
       read(line,*,iostat=ios) row
       if (ios >= 0) error stop 'Starting lattice has too many columns'
       spins(j+data%d,1+data%d:data%Nx+data%d) = int(row(1:data%Nx),ip)
    end do
    do
       read(u,'(a)',iostat=ios) line
       if (ios < 0) exit
       if (ios > 0) error stop 'Invalid lattice file'
       if (len_trim(line)>0) error stop 'Starting lattice has too many rows'
    end do
    close(u)
  end subroutine

  subroutine save_lattice(path,spins,data)
    character(len=*), intent(in) :: path
    type(basedata), intent(in) :: data
    integer(ip), intent(in) :: spins(:,:)
    integer :: u,j
    open(newunit=u,file=path,status='replace')
    do j=1,data%Ny
       write(u,*) spins(j+data%d,1+data%d:data%Nx+data%d)
    end do
    close(u)
  end subroutine

  subroutine shuffle_sites(sites)
    integer, intent(out) :: sites(:)
    integer :: i,j,tmp
    real(dp) :: r
    sites = [(i,i=1,size(sites))]
    do i=size(sites),2,-1
       call random_number(r)
       j=1+int(r*i)
       tmp=sites(i)
       sites(i)=sites(j)
       sites(j)=tmp
    end do
  end subroutine

  subroutine measure(spins,data,mean,interior,zero_fraction,maximum,interior_sites,map)
    integer(ip), intent(in) :: spins(:,:)
    type(basedata), intent(in) :: data
    real(dp), intent(out) :: mean,interior,zero_fraction,maximum
    integer, intent(out) :: interior_sites
    real(dp), intent(out), optional :: map(:,:)
    integer :: i,j,k,x,y,margin,zeros
    real(dp) :: e,norm
    mean=0
    interior=0
    maximum=0
    zeros=0
    interior_sites=0
    norm=sum(data%dval)
    ! Conservative full-neighborhood interior, using the actual nonzero stencil.
    margin=max(maxval(abs(data%xind-data%d-1)),maxval(abs(data%yind-data%d-1)))
    do i=1,data%Nx
       do j=1,data%Ny
          e=0
          do k=1,data%dnum
             x=i+data%xind(k)-1
             y=j+data%yind(k)-1
             if (spins(y,x)==spins(j+data%d,i+data%d)) e=e+data%dval(k)
          end do
          e=e/norm
          mean=mean+e
          maximum=max(maximum,e)
          ! Count exact zeros of the retained-weight model, not a chosen tolerance.
          if (e==0) zeros=zeros+1
          if (present(map)) map(j,i)=e
          if (i>margin .and. i<=data%Nx-margin .and. j>margin .and. j<=data%Ny-margin) then
             interior=interior+e
             interior_sites=interior_sites+1
          end if
       end do
    end do
    mean=mean/(real(data%Nx,dp)*data%Ny)
    zero_fraction=real(zeros,dp)/(real(data%Nx,dp)*data%Ny)
    if (interior_sites>0) then
       interior=interior/interior_sites
    else
       interior=ieee_value(0.0_dp,ieee_quiet_nan)
    end if
  end subroutine

  subroutine write_metadata(data)
    type(basedata), intent(in) :: data
    integer :: u
    open(newunit=u,file=trim(data%outfil)//'.meta',status='replace')
    write(u,'(a)') 'build_identity='//build_identity
    write(u,'(a)') 'compiler='//compiler_version()
    write(u,'(a)') 'compiler_options='//compiler_options()
    write(u,'(a,i0)') 'seed=',data%seed
    write(u,'(a,i0)') 'q=',data%q
    write(u,'(a,i0)') 'nx=',data%Nx
    write(u,'(a,i0)') 'ny=',data%Ny
    write(u,'(a,es25.17)') 'spacing=',data%a
    write(u,'(a)') 'kernel='//trim(data%kernel)
    write(u,'(a,l1)') 'random_order=',data%random_order
    write(u,'(a,i0)') 'completed_steps_before_run=',data%completed_steps
    write(u,'(a,i0)') 'temperature_steps=',data%tmax
    write(u,'(a,i0)') 'sweeps_before_offset=',data%ctmax
    write(u,'(a,i0)') 'sweep_offset=',data%offset
    write(u,'(a,es25.17)') 'cooling_factor=',data%tscale
    write(u,'(a,i0)') 'reheat_every=',data%reheat_every
    write(u,'(a,es25.17)') 'reheat_factor=',data%reheat_factor
    write(u,'(a)') 'initial_lattice='//trim(data%outfil)//'.initial.lat'
    write(u,'(a)') 'zero_definition=exact zero after discarding interaction weights <= 1e-14'
    close(u)
  end subroutine
end module
