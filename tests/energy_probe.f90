! Test driver: exercise production measurements and spin-energy differences.
program energy_probe
  use colorplane_consts
  use colorplane_utils
  use colorplane_reliability
  implicit none
  type(basedata) :: data
  integer(ip), allocatable :: spins(:,:)
  real(dp), allocatable :: map(:,:)
  real(dp) :: mean,inside,zero,maximum,after,eold,enew,weight,delta,t1,t2,t3
  integer :: n,i,j,k,x,y,cx,cy
  integer(ip) :: old,new
  character(len=1024) :: path
  if (updatenum(3,8,1)/=18) error stop 'Incorrect sweep count'
  read(*,*) data%q,data%Nx,data%Ny,data%a
  read(*,'(a)') data%kernel
  read(*,'(a)') path
  data%tmax=1
  data%ctmax=1
  data%tscale=0.9_dp
  call validate_input(data)
  call finalize_input(data)
  allocate(spins(data%Ny+2*data%d,data%Nx+2*data%d),map(data%Ny,data%Nx))
  spins=0
  call load_lattice(trim(path),spins,data)
  call measure(spins,data,mean,inside,zero,maximum,n,map)
  ! Compute the production-style local update delta at a central site.
  cx=data%d+data%Nx/2
  cy=data%d+data%Ny/2
  old=spins(cy,cx)
  new=newspin(old,data%q)
  eold=0
  enew=0
  do k=1,data%dnum
     x=cx+data%xind(k)-data%d-1
     y=cy+data%yind(k)-data%d-1
     weight=data%dval(k)
     if (spins(y,x)==old) eold=eold+weight
     if (spins(y,x)==new) enew=enew+weight
  end do
  spins(cy,cx)=new
  call measure(spins,data,after,t1,t2,t3,n)
  delta=2*(enew-eold) / &
       (real(data%Nx,dp)*data%Ny*sum(data%dval))
  ! Restore original and preserve the original map/measurements for Python.
  spins(cy,cx)=old
  call measure(spins,data,mean,inside,zero,maximum,n,map)
  print '(7es25.17)', mean,inside,zero,maximum,real(n,dp),after-mean,delta
  do j=1,data%Ny
     write(*,'(*(es25.17,1x))') (map(j,i),i=1,data%Nx)
  end do
end program
