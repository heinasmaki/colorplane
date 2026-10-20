! Independent-tile Metropolis sweeps on one core, OpenMP CPU, or OpenMP target.
! GPL-3.0-or-later.
module colorplane_parallel
  use colorplane_consts
  use colorplane_rng
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  !$ use omp_lib
  implicit none
  private
  public :: parallel_plan,make_plan,parallel_begin,parallel_reset,parallel_sweep
  public :: parallel_measure,parallel_finish,parallel_map,wall_seconds,configure_backend

  type :: parallel_plan
     integer :: tile=0,period=0,groups=0,stride=0,device=-1
     logical :: target_mode=.false.,threaded=.false.
     integer, allocatable :: starts(:),centers(:),sites(:),offsets(:),order(:)
     real(dp), allocatable :: weights(:),energy(:)
     integer(int64), allocatable :: accepted(:),changed(:)
  end type

contains
  function wall_seconds() result(t)
    real(dp) :: t
    integer(int64) :: ticks,rate
    call system_clock(ticks,rate)
    t=real(ticks,dp)/real(rate,dp)
  end function

  subroutine configure_backend(data,device)
    type(basedata), intent(inout) :: data
    integer, intent(out) :: device
    logical :: enabled,on_host
    integer :: devices
    enabled=.false.
    devices=0
    device=-1
    !$ enabled=.true.
    if (data%backend=='cpu' .or. data%backend=='gpu' .or. data%backend=='target-host') then
       if (.not. enabled) error stop 'This backend needs an OpenMP build; use build.py --openmp'
    end if
    if (data%threads>0) then
       if (data%backend/='cpu') error stop '--threads applies only to --backend cpu'
       !$ call omp_set_dynamic(.false.)
       !$ call omp_set_num_threads(data%threads)
    end if
    data%max_threads=1
    if (data%backend=='cpu') then
       !$ data%max_threads=omp_get_max_threads()
    end if
    if (data%backend=='gpu') then
       !$ devices=omp_get_num_devices()
       if (devices<1) error stop 'No OpenMP GPU device; refusing silent CPU fallback'
       !$ device=omp_get_default_device()
       if (device<0 .or. device>=devices) error stop 'Invalid OMP_DEFAULT_DEVICE'
       on_host=.true.
       !$omp target device(device) map(from:on_host)
       !$ on_host=omp_is_initial_device()
       !$omp end target
       if (on_host) error stop 'GPU offload fell back to host; GPU mode requires a real device'
    else if (data%backend=='target-host') then
       !$ device=omp_get_initial_device()
    end if
    data%device=device
  end subroutine

  subroutine make_plan(data,plan)
    type(basedata), intent(in) :: data
    type(parallel_plan), intent(out) :: plan
    integer, allocatable :: dx(:),dy(:),counts(:),next(:)
    integer :: radius,n,i,j,g,pos
    allocate(dx(data%dnum),dy(data%dnum))
    dx=data%xind-data%d-1
    dy=data%yind-data%d-1
    ! Within a tile all displacements have max(|dx|,|dy|) < tile.
    plan%tile=minval(max(abs(dx),abs(dy)))
    if (plan%tile<1) error stop 'Self-interaction invalidates independent-tile updates'
    radius=max(maxval(abs(dx)),maxval(abs(dy)))
    plan%period=(radius+plan%tile-1)/plan%tile+1
    ! Same-color tiles are separated by more than the stencil's bounding radius.
    if ((plan%period-1)*plan%tile+1<=radius) error stop 'Invalid tile separation'
    plan%groups=plan%period**2
    plan%stride=data%Ny+2*data%d
    plan%target_mode=data%backend=='gpu' .or. data%backend=='target-host'
    plan%threaded=data%backend=='cpu'
    plan%device=data%device
    n=data%Nx*data%Ny
    allocate(counts(plan%groups),next(plan%groups),plan%starts(plan%groups+1),plan%order(plan%groups))
    allocate(plan%centers(n),plan%sites(n),plan%offsets(data%dnum),plan%weights(data%dnum))
    allocate(plan%energy(n),plan%accepted(n),plan%changed(n))
    counts=0
    do i=1,data%Nx
       do j=1,data%Ny
          g=tile_group(i,j,plan%tile,plan%period)
          counts(g)=counts(g)+1
       end do
    end do
    plan%starts(1)=1
    do g=1,plan%groups
       plan%starts(g+1)=plan%starts(g)+counts(g)
    end do
    next=plan%starts(1:plan%groups)
    do i=1,data%Nx
       do j=1,data%Ny
          g=tile_group(i,j,plan%tile,plan%period)
          pos=next(g)
          plan%centers(pos)=j+data%d+(i+data%d-1)*plan%stride
          plan%sites(pos)=j+(i-1)*data%Ny
          next(g)=pos+1
       end do
    end do
    plan%offsets=dy+dx*plan%stride
    plan%weights=data%dval
    plan%order=[(g,g=1,plan%groups)]
    plan%accepted=0
    plan%changed=0
  end subroutine

  pure integer function tile_group(x,y,tile,period) result(group)
    integer, intent(in) :: x,y,tile,period
    group=1+mod((x-1)/tile,period)+period*mod((y-1)/tile,period)
  end function

  subroutine parallel_begin(s,plan)
    integer(ip), intent(inout) :: s(:)
    type(parallel_plan), intent(inout) :: plan
    if (plan%target_mode) call enter_device(s,plan%centers,plan%sites,plan%offsets,plan%weights, &
         plan%accepted,plan%changed,plan%energy,plan%device)
  end subroutine

  subroutine enter_device(s,centers,sites,offsets,weights,accepted,changed,energy,device)
    integer(ip), intent(inout) :: s(:)
    integer, intent(in) :: centers(:),sites(:),offsets(:),device
    real(dp), intent(in) :: weights(:)
    integer(int64), intent(inout) :: accepted(:),changed(:)
    real(dp), intent(inout) :: energy(:)
    !$omp target enter data device(device) map(to:s,centers,sites,offsets,weights) &
    !$omp& map(alloc:accepted,changed,energy)
  end subroutine

  subroutine parallel_reset(plan)
    type(parallel_plan), intent(inout) :: plan
    call reset_counts(plan%accepted,plan%changed,plan%target_mode,plan%threaded,plan%device)
  end subroutine

  subroutine reset_counts(accepted,changed,target_mode,threaded,device)
    integer(int64), intent(inout) :: accepted(:),changed(:)
    logical, intent(in) :: target_mode,threaded
    integer, intent(in) :: device
    integer :: p,n
    n=size(accepted)
    if (target_mode) then
       !$omp target teams distribute parallel do device(device) map(alloc:accepted,changed)
       do p=1,n
          accepted(p)=0
          changed(p)=0
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp parallel do if(threaded) schedule(static)
       do p=1,n
          accepted(p)=0
          changed(p)=0
       end do
       !$omp end parallel do
    end if
  end subroutine

  subroutine parallel_sweep(s,plan,data,temperature)
    integer(ip), intent(inout) :: s(:)
    type(parallel_plan), intent(inout) :: plan
    type(basedata), intent(in) :: data
    real(dp), intent(in) :: temperature
    integer :: c,g,first,last
    do c=1,plan%groups
       g=plan%order(c)
       first=plan%starts(g)
       last=plan%starts(g+1)-1
       if (last<first) cycle
       call update_group(s,plan%centers,plan%sites,plan%offsets,plan%weights,plan%accepted,plan%changed, &
            first,last,int(data%q),data%seed,data%completed_sweeps+1,temperature, &
            plan%target_mode,plan%threaded,plan%device)
    end do
  end subroutine

  subroutine update_group(s,centers,sites,offsets,weights,accepted,changed,first,last,q,seed,sweep,t, &
       target_mode,threaded,device)
    integer(ip), intent(inout) :: s(:)
    integer, intent(in) :: centers(:),sites(:),offsets(:),first,last,q,device
    real(dp), intent(in) :: weights(:),t
    integer(int64), intent(inout) :: accepted(:),changed(:)
    integer(int64), intent(in) :: seed,sweep
    logical, intent(in) :: target_mode,threaded
    integer :: p
    if (target_mode) then
       !$omp target teams distribute parallel do device(device) &
       !$omp& map(alloc:s,centers,sites,offsets,weights,accepted,changed)
       do p=first,last
          call advance_site(s,centers(p),sites(p),offsets,weights,q,seed,sweep,t,accepted(p),changed(p))
       end do
       !$omp end target teams distribute parallel do
    else
       !$omp parallel do if(threaded) schedule(static)
       do p=first,last
          call advance_site(s,centers(p),sites(p),offsets,weights,q,seed,sweep,t,accepted(p),changed(p))
       end do
       !$omp end parallel do
    end if
    ! Implicit barriers separate colors; interacting spins are never in flight together.
  end subroutine

  pure subroutine advance_site(s,center,site,offsets,weights,q,seed,sweep,t,accepted,changed)
    !$omp declare target
    integer(ip), intent(inout) :: s(:)
    integer, intent(in) :: center,site,offsets(:),q
    real(dp), intent(in) :: weights(:),t
    integer(int64), intent(in) :: seed,sweep
    integer(int64), intent(inout) :: accepted,changed
    integer :: candidate,k
    integer(ip) :: old,neighbor
    real(dp) :: eold,enew,u,delta
    logical :: take
    old=s(center)
    call random_move(seed,site,sweep,q,candidate,u)
    if (candidate>=old) candidate=candidate+1
    eold=0
    enew=0
    do k=1,size(offsets)
       neighbor=s(center+offsets(k))
       if (neighbor==old) eold=eold+weights(k)
       if (neighbor==candidate) enew=enew+weights(k)
    end do
    delta=enew-eold
    take=delta<=0
    if (.not. take) take=u<exp(-delta/t)
    if (take) then
       ! Same-kind assignment on NVIDIA's 32-bit build needs no INT(...,2)
       ! runtime helper on the device. GNU's compact build narrows 1..255 safely.
       s(center)=candidate
       accepted=accepted+1
       if (delta/=0) changed=changed+1
    end if
  end subroutine

  subroutine parallel_measure(s,plan,data,mean,inside,zero,maximum,ninside,accepted,changed)
    integer(ip), intent(in) :: s(:)
    type(parallel_plan), intent(inout) :: plan
    type(basedata), intent(in) :: data
    real(dp), intent(out) :: mean,inside,zero,maximum
    integer, intent(out) :: ninside
    integer(int64), intent(out) :: accepted,changed
    integer :: site,i,j,margin,zeros
    real(dp) :: e
    call energy_kernel(s,plan%offsets,plan%weights,plan%energy,plan%accepted,plan%changed, &
         data%Nx,data%Ny,data%d,plan%stride,plan%target_mode,plan%threaded,plan%device)
    mean=0
    inside=0
    maximum=0
    zeros=0
    ninside=0
    margin=max(maxval(abs(data%xind-data%d-1)),maxval(abs(data%yind-data%d-1)))
    ! Fixed host summation order makes CPU diagnostics independent of thread count.
    do site=1,size(plan%energy)
       i=1+(site-1)/data%Ny
       j=1+mod(site-1,data%Ny)
       e=plan%energy(site)
       mean=mean+e
       maximum=max(maximum,e)
       if (e==0) zeros=zeros+1
       if (i>margin .and. i<=data%Nx-margin .and. j>margin .and. j<=data%Ny-margin) then
          inside=inside+e
          ninside=ninside+1
       end if
    end do
    mean=mean/size(plan%energy)
    zero=real(zeros,dp)/size(plan%energy)
    if (ninside>0) then
       inside=inside/ninside
    else
       inside=ieee_value(0.0_dp,ieee_quiet_nan)
    end if
    accepted=sum(plan%accepted)
    changed=sum(plan%changed)
  end subroutine

  subroutine energy_kernel(s,offsets,weights,energy,accepted,changed,nx,ny,d,stride,target_mode,threaded,device)
    integer(ip), intent(in) :: s(:)
    integer, intent(in) :: offsets(:),nx,ny,d,stride,device
    real(dp), intent(in) :: weights(:)
    real(dp), intent(inout) :: energy(:)
    integer(int64), intent(inout) :: accepted(:),changed(:)
    logical, intent(in) :: target_mode,threaded
    integer :: site,center,k,n
    real(dp) :: norm,e
    n=nx*ny
    norm=sum(weights)
    if (target_mode) then
       !$omp target teams distribute parallel do device(device) &
       !$omp& map(alloc:s,offsets,weights,energy) private(center,k,e)
       do site=1,n
          center=1+mod(site-1,ny)+d+((site-1)/ny+d)*stride
          e=0
          do k=1,size(offsets)
             if (s(center+offsets(k))==s(center)) e=e+weights(k)
          end do
          energy(site)=e/norm
       end do
       !$omp end target teams distribute parallel do
       !$omp target update device(device) from(energy,accepted,changed)
    else
       !$omp parallel do if(threaded) schedule(static) private(center,k,e)
       do site=1,n
          center=1+mod(site-1,ny)+d+((site-1)/ny+d)*stride
          e=0
          do k=1,size(offsets)
             if (s(center+offsets(k))==s(center)) e=e+weights(k)
          end do
          energy(site)=e/norm
       end do
       !$omp end parallel do
    end if
  end subroutine

  subroutine parallel_finish(s,plan)
    integer(ip), intent(inout) :: s(:)
    type(parallel_plan), intent(inout) :: plan
    if (plan%target_mode) call leave_device(s,plan%centers,plan%sites,plan%offsets,plan%weights, &
         plan%accepted,plan%changed,plan%energy,plan%device)
  end subroutine

  subroutine leave_device(s,centers,sites,offsets,weights,accepted,changed,energy,device)
    integer(ip), intent(inout) :: s(:)
    integer, intent(in) :: centers(:),sites(:),offsets(:),device
    real(dp), intent(in) :: weights(:)
    integer(int64), intent(inout) :: accepted(:),changed(:)
    real(dp), intent(inout) :: energy(:)
    !$omp target exit data device(device) map(from:s) &
    !$omp& map(delete:centers,sites,offsets,weights,accepted,changed,energy)
  end subroutine

  subroutine parallel_map(plan,data)
    type(parallel_plan), intent(in) :: plan
    type(basedata), intent(in) :: data
    integer :: u,j,i
    open(newunit=u,file=trim(data%outfil)//'.ngb',status='replace')
    do j=1,data%Ny
       write(u,*) (plan%energy(j+(i-1)*data%Ny),i=1,data%Nx)
    end do
    close(u)
  end subroutine
end module
