subroutine adjust_omega_BBH(flag_omega,error_ome)
  use grid_parameter, only : mass_eps, eps
  use def_quantities, only : admmass_asymp, komarmass_asymp
  use def_bh_parameter, only : ome_bh, spin_bh
  implicit none
  integer :: flag_omega
  integer, save :: count_adj, ii
  real(8), save :: fnc_val, fnc_val_bak, dfdval, &
  &                ome_bh_bak, delta_ome_bh, delta_ome_bh_bak
  real(8), save :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
!  real(8), save :: facp(5) = (/ 1.0d-1, 2.0d-1, 3.0d-1, 4.0d-1, 5.0d-1 /)
  real(8) :: error_ome, fac0 = 1.0d0
!  real(8), save :: facp(10) = (/ 2.0d-2, 4.0d-2, 6.0d-2, 8.0d-2, &
!  &              1.0d-1, 2.0d-1, 4.0d-1, 6.0d-1, 8.0d-1, 1.0d-0 /)
!
! -- secant method: iteration for omega for a circular solution
!
  if (flag_omega.eq.99) then 
    if (admmass_asymp.le.komarmass_asymp) delta_ome_bh = -0.05d0*ome_bh
    if (admmass_asymp.ge.komarmass_asymp) delta_ome_bh =  0.05d0*ome_bh
    ome_bh = ome_bh + delta_ome_bh
    fnc_val = (admmass_asymp - komarmass_asymp)/admmass_asymp
    flag_omega = 1
    count_adj = 0
    return
  end if
  fnc_val_bak = fnc_val
  delta_ome_bh_bak = delta_ome_bh
  fnc_val = (admmass_asymp - komarmass_asymp)/admmass_asymp
  dfdval = (fnc_val - fnc_val_bak)/delta_ome_bh
  delta_ome_bh = - fnc_val/dfdval
!
  error_ome = dabs(delta_ome_bh/ome_bh)
!
  if (error_ome.le.eps) then 
    fnc_val =  fnc_val_bak
    delta_ome_bh = delta_ome_bh_bak
    flag_omega = 0
    return
  end if
!
  count_adj = count_adj + 1
  ii = min0(5,count_adj)
  fac0 = facp(ii)
  delta_ome_bh = fac0*delta_ome_bh
  ome_bh_bak = ome_bh
  ome_bh = ome_bh + delta_ome_bh
!
    write(6,'(1a24,2(7x,i5))') ' -- ome iteration # --  ',count_adj
    write(6,'(1a24,1p,2e12.4)')' -- ome OLD -> NEW  --  ',ome_bh_bak, ome_bh
    write(6,*)'-- ADM and Komar --',  admmass_asymp,  komarmass_asymp
    write(6,*)'-- 1 - MK/MADM   --',  fnc_val_bak,  fnc_val
!
!!write(6,*)flag_restmass, count_adj
!!write(6,*)ii, fac0
!!write(6,*)fixeddlm,fnc_val,hi
!!write(6,*)restmass, restmass_sph
!!write(6,*)admmass,  komarmass
!!write(6,*)emdc,par_val
!!write(6,*)sgg,gg
!!write(6,*)epsmax_fnc, hi_fac
!
  return
end subroutine adjust_omega_BBH
