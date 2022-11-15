SUBROUTINE quark_h2rho(quark_rho2h,quark_rho2h_dot,rho1,rho2,quarkh,rho)
  IMPLICIT NONE
  interface
      function quark_rho2h(rho)
      real(8), intent(in) :: rho
!      integer :: ii
      real(8) :: quark_rho2h
      end function quark_rho2h
      function quark_rho2h_dot(rho)
      real(8), intent(in) :: rho
!      integer :: ii
      real(8) :: quark_rho2h_dot
      end function
  end interface
  real(8), intent(in) :: rho1, rho2, quarkh
  real(8), intent(inout) :: rho
  real(8) :: quarkh_tmp, rhomax, rhomin
  integer :: time
!  real(8), external :: quark_rho2p
!  real(8), save :: small
!  real(8) :: pgamma
!  open(3,file='test1.dat')
  time = 0
!  pgamma = 13./3.
!  if (quarkh.le.1.0d0) then
!     rho = ((pgamma-1.0d0)/pgamma)*(quarkh-quark_rho2h(small))+ &
!          & quark_rho2p(small)/small
!     if (rho.gt.0.0d0) small = rho
!     return
!  end if
!------------------------------------------------------------------------
! These commented lines are written when I wanted to test my code with NS
!EoS states. In practice these will never be needed any more if the _qeos
! codes are only used for quark stars. But in case anyone has made a really 
! bad decision of using this _qeos codes to build neutron stars, then you 
!will need to uncomment the lines above and put pgamma to be the gamma index
! of the peos that you use, and also you need to uncomment the lines in 
!../../Subroutine/calc_surface_qeos.f90.
!                         Enping Zhou  July 2016

!------------------------------------------------------------------------
  quarkh_tmp = quark_rho2h(rho)
  do
!    quarkh_tmp = quark_rho2h(rho)
    rho = rho + (quarkh-quarkh_tmp)/quark_rho2h_dot(rho)
    quarkh_tmp = quark_rho2h(rho)
!    write(3,*) quarkh_tmp, quarkh,rho
    if (abs(quarkh_tmp-quarkh)/quarkh.lt.1e-10) then
!       write(3,*) time+1
!     if (rho.gt.0.0d0) small = rho
       return 
    end if
    time = time +1
    if (time.gt.20) then
!       write(3,*) 'newton method failed'
       exit
    end if
  end do
  time = 0
  rhomin = rho1
  rhomax = rho2
  do 
    rho = 0.5d0*(rhomin+rhomax)
    quarkh_tmp =  quark_rho2h(rho)
    if (abs(quarkh_tmp-quarkh)/quarkh.lt.1e-10) then
!       if (rho.gt.0.0d0) small = rho
       return
    end if
    if (quarkh_tmp.gt.quarkh) then
       rhomax = rho
    end if
    if (quarkh_tmp.lt.quarkh) then
       rhomin = rho
    end if 
!    time = time + 1 
!    write(3,*) time, quark_rho2h(rho), quarkh, quarkh1, quarkh2, rhomin, rhomax
  end do
end SUBROUTINE quark_h2rho
