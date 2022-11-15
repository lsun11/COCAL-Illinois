subroutine adjust_charge(flag_st,charge_type)
  use def_quantities, only : charge, charge_asymp
  use integrability_fnc_MHD, only : MHDpar_charge
  implicit none
  integer :: flag_st, ii
  integer, save :: count_st
  real(8), save :: x_st(0:2), f_st(0:1)
  real(8) :: delta, dfdx, error_st
  real(8) :: charge_adjust, charge_hope = 0.0d0, small = 1.0d-10 
  real(8) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
  character(len=6) :: charge_type
!
! -- secant method: iteration for charge to be zero
!
  if (charge_type.eq.'asympt') then
    call calc_charge_asympto('ns')
    charge_adjust = charge_asymp
  else if (charge_type.eq.'volume') then
    call calc_charge_MHD
    charge_adjust = charge
  else
    write(6,*) 'charge type wrong'
    stop
  end if
!
  if (flag_st.ne.0) then 
    x_st(0) = x_st(1)
    f_st(0) = f_st(1)
  end if
!
  error_st = charge_adjust - charge_hope
  write(6,'(1a24,2(7x,i5))') ' -- charge flag   #  -- ',flag_st
  write(6,'(1a24,1p,2e12.4)')' -- charge and param -- ',charge_adjust, &
  &                                                     MHDpar_charge
!
  if (dabs(error_st).gt.small) then 
    x_st(1) = MHDpar_charge
    f_st(1) = error_st
  else 
    write(6,'(1a24)') ' -- charge converged -- '
    return 
  end if
!
  if (flag_st.eq.0) then 
    count_st = 1
    delta = MHDpar_charge
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
  MHDpar_charge = x_st(2)
!
  if (flag_st.eq.0) flag_st = +1
!
    write(6,'(1a6,1p,3e12.4)')'charge', x_st(0),x_st(1),x_st(2)
    write(6,'(1a6,1p,3e12.4)')'charge', f_st(0),f_st(1)
    write(6,'(1a6,1p,3e12.4)')'charge', delta, dfdx
!    write(6,*)'-- ADM and Komar --',  admmass_asymp,  komarmass_asymp
!    write(6,*)'-- 1 - MK/MADM   --',  fnc_val_bak,  fnc_val
!
  return
end subroutine adjust_charge
