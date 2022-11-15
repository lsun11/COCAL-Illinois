subroutine search_emdmax_xaxis(iremax,emdmax,xsol)
  use phys_constant,  only  : long
  use grid_parameter, only  : nrf, ntfxy
  use coordinate_grav_r, only  : rg
  use def_matter, only  :   emd
  implicit none
  real(long) :: emdmax, xsol, minB, niA, sqrtD
  real(long) :: emdmaxp, emdmaxm, xsolp, xsolm
  real(long) :: x(4), f(4), x1, x2, x3, x4, f1, f2, f3, f4, &
  &             x12, x13, x14, x21, x23, x24, &
  &             x31, x32, x34, x41, x42, x43
  integer    :: ir, iremax, ir0
  real(long), external :: lagint_4th
!
! search maximum emd on positive x-axis
  if (iremax.eq.0) then
    emdmax = emd(iremax,ntfxy,0)
    xsol = rg(0)
    return
  else
    ir0 = iremax-1
    if (emd(iremax-1,ntfxy,0).gt.emd(iremax+1,ntfxy,0)) ir0=max0(0,iremax-2)
    x(1:4) =  rg(ir0:ir0+3)
    f(1:4) = emd(ir0:ir0+3,ntfxy,0)
    x1 = x(1) ; x2 = x(2) ; x3 = x(3) ; x4 = x(4)
    f1 = f(1) ; f2 = f(2) ; f3 = f(3) ; f4 = f(4)
    x12 = x(1)-x(2)
    x13 = x(1)-x(3)
    x14 = x(1)-x(4)
    x21 = x(2)-x(1)
    x23 = x(2)-x(3)
    x24 = x(2)-x(4)
    x31 = x(3)-x(1)
    x32 = x(3)-x(2)
    x34 = x(3)-x(4)
    x41 = x(4)-x(1)
    x42 = x(4)-x(2)
    x43 = x(4)-x(3)
!
minB =  f4*x1*x12*x13*x14*x21*x23*x24*x31*x32*x34 + &
     &  f4*x12*x13*x14*x2*x21*x23*x24*x31*x32*x34 + &
     &  f4*x12*x13*x14*x21*x23*x24*x3*x31*x32*x34 + &
     &  f3*x1*x12*x13*x14*x21*x23*x24*x41*x42*x43 + &
     &  f3*x12*x13*x14*x2*x21*x23*x24*x41*x42*x43 + &
     &  f2*x1*x12*x13*x14*x31*x32*x34*x41*x42*x43 + &
     &  f1*x2*x21*x23*x24*x31*x32*x34*x41*x42*x43 + &
     &  f2*x12*x13*x14*x3*x31*x32*x34*x41*x42*x43 + &
     &  f1*x21*x23*x24*x3*x31*x32*x34*x41*x42*x43 + &
     &  f3*x12*x13*x14*x21*x23*x24*x4*x41*x42*x43 + &
     &  f2*x12*x13*x14*x31*x32*x34*x4*x41*x42*x43 + &
     &  f1*x21*x23*x24*x31*x32*x34*x4*x41*x42*x43 
!
sqrtD=dsqrt((f4*x12*x13*x14*x21*x23*x24*(x1 + x2 + x3)*x31*x32*x34 + &
     &       (f3*x12*x13*x14*x21*x23*x24*(x1 + x2 + x4) + &
     &          x31*x32*x34*(f2*x12*x13*x14*(x1 + x3 + x4) + &
     &             f1*x21*x23*x24*(x2 + x3 + x4)))*x41*x42*x43)**2 - &
     &    3*(f4*x12*x13*x14*x21*x23*x24*x31*x32*x34 + &
     &       (f3*x12*x13*x14*x21*x23*x24 + &
     &          (f2*x12*x13*x14 + f1*x21*x23*x24)*x31*x32*x34)*x41*x42*x43) &
     &      *(f4*x12*x13*x14*x21*x23*x24*(x2*x3 + x1*(x2 + x3))*x31*x32* &
     &        x34 + (f3*x12*x13*x14*x21*x23*x24*(x2*x4 + x1*(x2 + x4)) + &
     &          x31*x32*x34*(f2*x12*x13*x14*(x3*x4 + x1*(x3 + x4)) + &
     &             f1*x21*x23*x24*(x3*x4 + x2*(x3 + x4))))*x41*x42*x43))
!
niA  = 3.*(f4*x12*x13*x14*x21*x23*x24*x31*x32*x34 + &
     &    (f3*x12*x13*x14*x21*x23*x24 + &
     &       (f2*x12*x13*x14 + f1*x21*x23*x24)*x31*x32*x34)*x41*x42*x43)
!
    xsolp = (minB + sqrtD)/niA
    xsolm = (minB - sqrtD)/niA
!
    emdmaxp = lagint_4th(x,f,xsolp)
    emdmaxm = lagint_4th(x,f,xsolm)
!
    emdmax = emdmaxp
    xsol   = xsolp
    if (emdmaxm.gt.emdmaxp) then 
      emdmax = emdmaxm
      xsol   = xsolm
    end if
!
  end if
!
end subroutine search_emdmax_xaxis
