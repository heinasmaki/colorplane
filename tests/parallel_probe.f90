! Test driver for the coloring, counter RNG, and complete blocked sweeps.
program parallel_probe
  use colorplane_consts
  use colorplane_rng
  use colorplane_parallel
  use colorplane_utils
  use colorplane_reliability
  use, intrinsic :: iso_fortran_env, only: int64
  implicit none
  type(basedata) :: data
  type(parallel_plan) :: plan
  integer(ip), allocatable, target :: spins(:,:)
  integer(ip), pointer, contiguous :: flat(:)
  integer(int64) :: counter(4),key(2),words(4),seed,sweep,accepted,changed
  integer :: ios,site,q,candidate,i,j,g,p,n,device,loops
  real(dp) :: u,t,mean,inside,zero,maximum
  character(len=32) :: mode
  call get_command_argument(1,mode)
  select case(trim(mode))
  case('rng')
     do
        read(*,*,iostat=ios) counter,key
        if (ios<0) exit
        if (ios/=0) error stop 'Invalid RNG probe input'
        call philox4x32(counter,key,words)
        print '(4(z8.8,1x))',words
     end do
  case('moves')
     do
        read(*,*,iostat=ios) seed,site,sweep,q
        if (ios<0) exit
        if (ios/=0) error stop 'Invalid move probe input'
        call random_move(seed,site,sweep,q,candidate,u)
        print '(i0,1x,es25.17)',candidate,u
     end do
  case('plan','sweep')
     read(*,*) data%q,data%Nx,data%Ny,data%a
     read(*,'(a)') data%kernel
     data%tmax=1
     data%ctmax=1
     data%tscale=0.9_dp
     call validate_input(data)
     call finalize_input(data)
     if (mode=='plan') then
        call make_plan(data,plan)
        print '(4(i0,1x))',plan%tile,plan%period,plan%groups,data%d
        do g=1,plan%groups
           do p=plan%starts(g),plan%starts(g+1)-1
              print '(3(i0,1x))',plan%sites(p),g,plan%centers(p)
           end do
        end do
     else
        call get_command_argument(2,data%backend)
        read(*,*) data%seed,t,loops
        call configure_backend(data,device)
        call make_plan(data,plan)
        allocate(spins(data%Ny+2*data%d,data%Nx+2*data%d))
        spins=0
        do i=1,data%Nx
           do j=1,data%Ny
              spins(j+data%d,i+data%d)=int(1+mod(3*i+7*j+i*j,int(data%q)),ip)
           end do
        end do
        flat(1:size(spins))=>spins
        call parallel_begin(flat,plan)
        call parallel_reset(plan)
        do i=1,loops
           call parallel_sweep(flat,plan,data,t)
           data%completed_sweeps=data%completed_sweeps+1
        end do
        call parallel_measure(flat,plan,data,mean,inside,zero,maximum,n,accepted,changed)
        call parallel_finish(flat,plan)
        print '(es25.17,1x,2(i0,1x))',mean,accepted,changed
        do j=1,data%Ny
           write(*,'(*(i0,1x))') spins(j+data%d,1+data%d:data%Nx+data%d)
        end do
     end if
  case default
     error stop 'Unknown probe mode'
  end select
end program
