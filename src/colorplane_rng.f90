! Counter-based random numbers for independent, schedule-invariant spin updates.
! GPL-3.0-or-later. Independent Fortran implementation of Philox4x32-10.
! Algorithm/constants: D. E. Shaw Research, Random123; see THIRD_PARTY_NOTICES.md.
module colorplane_rng
  use, intrinsic :: iso_fortran_env, only: int64
  use colorplane_consts, only: dp
  implicit none
  private
  public :: philox4x32, random_move
  integer(int64), parameter :: mask32=int(z'FFFFFFFF',int64), mask16=65535_int64
contains
  pure subroutine multiply32(a,b,hi,lo)
    !$omp declare target
    integer(int64), intent(in) :: a,b
    integer(int64), intent(out) :: hi,lo
    integer(int64) :: lower,upper
    ! Two 16x32 products avoid undefined signed 64-bit integer overflow.
    lower=iand(a,mask16)*b
    upper=shiftr(a,16)*b+shiftr(lower,16)
    lo=ior(shiftl(iand(upper,mask16),16),iand(lower,mask16))
    hi=shiftr(upper,16)
  end subroutine

  pure subroutine philox4x32(counter,key,result)
    !$omp declare target
    integer(int64), intent(in) :: counter(4),key(2)
    integer(int64), intent(out) :: result(4)
    integer(int64) :: c(4),k(2),h0,l0,h1,l1
    integer :: round
    c=counter
    k=key
    do round=1,10
       call multiply32(3528531795_int64,c(1),h0,l0)
       call multiply32(3449720151_int64,c(3),h1,l1)
       c=[ieor(ieor(h1,c(2)),k(1)),l1,ieor(ieor(h0,c(4)),k(2)),l0]
       k(1)=iand(k(1)+2654435769_int64,mask32)
       k(2)=iand(k(2)+3144134277_int64,mask32)
    end do
    result=c
  end subroutine

  pure subroutine random_move(seed,site,sweep,q,candidate,u)
    !$omp declare target
    integer(int64), intent(in) :: seed,sweep
    integer, intent(in) :: site,q
    integer, intent(out) :: candidate
    real(dp), intent(out) :: u
    integer(int64) :: counter(4),key(2),words(4),limit,n
    n=int(q-1,int64)
    limit=(4294967296_int64/n)*n
    key=[iand(seed,mask32),shiftr(seed,32)]
    counter=[int(site,int64),iand(sweep,mask32),shiftr(sweep,32),0_int64]
    do
       call philox4x32(counter,key,words)
       if (words(1)<limit) exit
       counter(4)=counter(4)+1
    end do
    ! Rejection sampling gives an unbiased choice among q-1 alternative colors.
    candidate=1+int(modulo(words(1),n))
    ! Exactly representable 53-bit uniform in [0,1), independent of proposal word.
    u=(real(words(3),dp)*2097152.0_dp+real(shiftr(words(4),11),dp))*2.0_dp**(-53)
  end subroutine
end module
