subroutine calc_eqsolver_TrpBH(rr,rini,rval)
  use phys_constant, only : long
  use def_bh_parameter, only : mass_pBH
  implicit none
  integer :: count, ii
  real(long) :: rr, rini, rval
  real(long) :: R1, R2, F1, F2, Rtmp, M, jacobian, fac0, error
  real(long) :: facp(5) = (/ 2.0d-1, 5.0d-1, 1.0d-0, 1.0d-0, 1.0d-0 /)
!
! --- Solve parametric trunpet BH eqs using secant method.
!
  M = mass_pBH
  R1 = rini
  R2 = rini*1.05d0
  F1 = fnc_R_trpBH(R1,M) - rr
!
  count = 0
  do 
    count = count + 1
    ii = min0(5,count)
    fac0 = facp(ii)
!
    F2 = fnc_R_trpBH(R2,M) - rr
    jacobian = (F2-F1)/(R2-R1)
    Rtmp = R2 - F2/jacobian
    R1 = R2
    F1 = F2
!
    R2 = fac0*Rtmp + (1.0d0-fac0)*R2
    if (R2.le.3.0d0*M/2.0d0) R2 = 3.0d0*M/2.0d0
    error = 2.d0*dabs(R2 - R1)/(dabs(R2) + dabs(R1))
! 
    IF (error <= 1.d-14) EXIT 
    IF (count.ge.1000) exit
  end do
  IF (count.ge.1000) stop 'TrpBH failed'
!
  rval   = R2
!
  contains
    function fnc_R_trpBH(R,M)
      real(long) :: R, M, fnc_R_trpBH
      fnc_R_trpBH = 0.25d0*(2.0d0*R+M+dsqrt(4.0d0*R**2+4.0d0*M*R+3.0d0*M**2)) &
      &           *((4.0d0+3.0d0*dsqrt(2.0d0))*(2.0d0*R-3.0d0*M) &
      &  /(8.0d0*R+6.0d0*M+3.0d0*dsqrt(8.0d0*R**2+8.0d0*M*R+6.0d0*M**2)) &
      &            )**(1.0d0/dsqrt(2.0d0))
    end function fnc_R_trpBH
end subroutine calc_eqsolver_TrpBH
