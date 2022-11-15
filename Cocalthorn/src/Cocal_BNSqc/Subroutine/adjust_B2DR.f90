subroutine adjust_B2DR(flag_st)
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntfeq, npfxzp
  use def_matter, only : omef
  use def_matter_parameter, only : ome, B2DR, DRAT_B2DR
  use make_array_1d
  use interface_search_max_radial
  implicit none
  integer :: flag_st, ii
  integer, save :: count_st
  real(long), save :: x_st(0:2), f_st(0:1)
  real(long) :: delta, dfdx, error_st, diff_st
  real(long) :: DRAT_adjust, DRAT_hope, small = 1.0d-10 
  real(long) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
  character(len=6) :: charge_type
  real(long), pointer :: omef_x(:)
  real(long)          :: omef_max, x_omax
  integer             :: irf_omax
!
! -- secant method: iteration for charge to be zero
!
  call alloc_array1d(omef_x,0,nrf)
  omef_x(0:nrf)=omef(0:nrf,ntfeq,npfxzp)
  irf_omax = maxloc(omef_x,DIM=1) - 1
!  write(6,*) maxloc(omef_x) - 1
!  write(6,*)'ome max',omef_x(irf_omax-1),omef_x(irf_omax),omef_x(irf_omax+1)
!stop
  call search_max_radial(irf_omax,omef_x,omef_max,x_omax)
  deallocate(omef_x)
  DRAT_adjust = omef_max/ome
  DRAT_hope   = DRAT_B2DR
!
  if (flag_st.ne.0) then 
    x_st(0) = x_st(1)
    f_st(0) = f_st(1)
  end if
!
  diff_st = DRAT_adjust - DRAT_hope
  error_st = dabs(diff_st/DRAT_hope)
  write(6,'(1a20,2(7x,i5))') '-- DRAT flag   # -- ',flag_st
  write(6,'(1a20,1p,2e12.4)')'-- DRAT and B2DR -- ',DRAT_adjust, B2DR
!
  if (error_st.gt.small) then 
    x_st(1) = B2DR
    f_st(1) = diff_st
  else 
    write(6,'(1a24)') ' -- B2DR converged -- '
    return 
  end if
!
  if (flag_st.eq.0) then 
    count_st = 1
    delta = 0.1*B2DR
    if (delta.eq.0.0d0) delta = 0.1d0
  else
    count_st = count_st + 1
    dfdx = (f_st(1) - f_st(0))/(x_st(1) - x_st(0))
    delta = - f_st(1)/dfdx
    flag_st = +1
    if (delta.lt.0.0d0) flag_st = -1
  end if
!
  ii = min0(5,count_st)
  x_st(2) = x_st(1) + delta*facp(ii)
  B2DR = x_st(2)
!
  if (flag_st.eq.0) flag_st = +1
!
    write(6,'(1a6,1p,3e12.4)')'B2DR', x_st(0),x_st(1),x_st(2)
    write(6,'(1a6,1p,3e12.4)')'B2DR', f_st(0),f_st(1)
    write(6,'(1a6,1p,3e12.4)')'B2DR', delta, dfdx
!
  return
end subroutine adjust_B2DR
