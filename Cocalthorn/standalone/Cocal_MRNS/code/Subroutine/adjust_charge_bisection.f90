subroutine adjust_charge_bisection(flag_bis,charge_type)
  use grid_parameter, only : mass_eps, eps
  use def_quantities, only : charge, charge_asymp
  use integrability_fnc_MHD, only : MHDpar_charge
  implicit none
  integer :: flag_bis
  real(8), save :: delta, x_bis(0:1), f_bis(0:1)
  real(8) :: error_bis, charge_adjust, charge_hope = 0.0d0, small = 1.0d-10 
  character(len=6) :: charge_type
!
! -- bisection method: iteration for charge to be zero
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
  if (flag_bis.ne.0) then 
    x_bis(0) = x_bis(1)
    f_bis(0) = f_bis(1)
  end if
!
  error_bis = charge_adjust - charge_hope
  write(6,'(1a24,2(7x,i5))') ' -- charge flag   #  -- ',flag_bis
  write(6,'(1a24,1p,2e12.4)')' -- charge and param -- ',charge_adjust, &
  &                                                     MHDpar_charge
!
  if (dabs(error_bis).gt.small) then 
    x_bis(1) = MHDpar_charge
    f_bis(1) = error_bis
  else 
    write(6,'(1a24)') ' -- charge converged -- '
    return 
  end if
!
  if (flag_bis.eq.0) then 
    delta = MHDpar_charge*0.5d0
    if (delta.eq.0.0d0) delta = 0.1d0
    flag_bis = +1
  end if
!
  if (flag_bis.ne.0) then
    if (f_bis(1)*f_bis(0).ge.0) then 
      if (dabs(f_bis(1)).gt.dabs(f_bis(0))) then 
        delta    = - delta
        flag_bis = - flag_bis
      end if
    else 
        delta    = - delta
        flag_bis = - flag_bis
        if (delta.gt.0.0d0.and.dabs(error_bis).gt.small) delta = 0.5d0*delta
    end if
  end if
!
  MHDpar_charge = MHDpar_charge + delta
!
    write(6,'(1a6,1p,3e12.4)')'charge', x_bis(0),x_bis(1)
    write(6,'(1a6,1p,3e12.4)')'charge', f_bis(0),f_bis(1)
    write(6,'(1a6,1p,3e12.4)')'charge', delta
!    write(6,*)'-- ADM and Komar --',  admmass_asymp,  komarmass_asymp
!    write(6,*)'-- 1 - MK/MADM   --',  fnc_val_bak,  fnc_val
!
  return
end subroutine adjust_charge_bisection
