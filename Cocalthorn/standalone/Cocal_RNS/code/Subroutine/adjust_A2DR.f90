subroutine adjust_A2DR(flag_st,istep_st)
  use grid_parameter, only : nrf, ntfeq, npfxzp, eps
  use def_matter, only : omef
  use def_matter_parameter, only : ome, A2DR, DRAT_A2DR
  implicit none
  integer :: flag_st, ii, istep_st
  integer, save :: count_st
  real(8), save :: x_st(0:2), f_st(0:1)
  real(8) :: delta, dfdx, error_st, diff_st
  real(8) :: DRAT_adjust, DRAT_hope, small = 1.0d-8
  real(8) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
  character(len=6) :: charge_type
!
! -- secant method: iteration for charge to be zero
!
  DRAT_adjust = omef(nrf,ntfeq,npfxzp)/ome
  DRAT_hope   = DRAT_A2DR
!
  if (flag_st.ne.0) then 
    x_st(0) = x_st(1)
    f_st(0) = f_st(1)
  end if
!
  diff_st = DRAT_adjust - DRAT_hope
  error_st = dabs(diff_st/DRAT_hope)
  write(6,'(1a20,2(7x,i5))') '-- DRAT flag   # -- ',flag_st
  write(6,'(1a20,1p,2e12.4)')'-- DRAT and A2DR -- ',DRAT_adjust, A2DR
!
  if (error_st.gt.eps) then 
    x_st(1) = A2DR
    f_st(1) = diff_st
  else 
    write(6,'(1a24)') ' -- A2DR converged -- '
    flag_st = 0
    return 
  end if
!
  if (istep_st.eq.-1) then 
    istep_st = 0
    count_st = 1
    delta = 0.1*A2DR
    if (delta.eq.0.0d0) delta = 0.1d0
  else
    count_st = count_st + 1
    dfdx = (f_st(1) - f_st(0))/(x_st(1) - x_st(0))
    delta = - f_st(1)/dfdx
  end if
!
  ii = min0(5,count_st)
  x_st(2) = x_st(1) + delta*facp(ii)
  A2DR = x_st(2)
!
  flag_st = +1
!
    write(6,'(1a6,1p,3e12.4)')'A2DR', x_st(0),x_st(1),x_st(2)
    write(6,'(1a6,1p,3e12.4)')'A2DR', f_st(0),f_st(1)
    write(6,'(1a6,1p,3e12.4)')'A2DR', delta, dfdx
!
  return
end subroutine adjust_A2DR
