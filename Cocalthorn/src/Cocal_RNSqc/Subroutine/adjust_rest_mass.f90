subroutine adjust_rest_mass(flag_restmass,count_adj)
  use grid_parameter, only : mass_eps
  use def_matter_parameter, only : emdc
  use def_quantities
  implicit none
  integer :: flag_restmass
  integer :: count_adj, ii
  real(8), save :: hi, hi_fac, par_val, fnc_val, sgg, gg, &
  &                epsmax_fnc, epsmax_par, fac0
  real(8), save :: fixeddlm, fixedvir
  real(8), save :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
!  real(8), save :: facp(10) = (/ 2.0d-2, 4.0d-2, 6.0d-2, 8.0d-2, &
!  &              1.0d-1, 2.0d-1, 4.0d-1, 6.0d-1, 8.0d-1, 1.0d-0 /)
!
! -- Discrete Newton method: iteration for a constant rest mass.
!
  fixeddlm = (restmass - restmass_sph)/restmass_sph
  fixedvir = (admmass  -    komarmass)/admmass
  epsmax_fnc = dabs(fixeddlm)
!
  if(flag_restmass.eq.0) then 
    epsmax_par = 1.0d0
    if (count_adj.ne.0) epsmax_par = dabs((emdc - par_val)/emdc)
    par_val = emdc
    fnc_val = restmass
!    fnc_val = fixeddlm 
    hi_fac = - dmin1(1.0d-02,10.0d0**dble(int(dlog10(epsmax_par))-1))
!!   hi_fac = dmin1(1.0d-02,10.0d0**dble(int(dlog10(epsmax_par))-1))
!    hi_fac = dmin1(1.0d-02,10.0d0**dble(int(dlog10(epsmax_par))))
    hi = hi_fac*dabs(par_val)
    emdc = par_val + hi
!wrong?    emdc = par_val - hi
    flag_restmass = 1
    return
  end if
  if(flag_restmass.ge.1) then 
    count_adj = count_adj + 1 
!    sgg= - fnc_val
!    gg = (fixeddlm - fnc_val)/hi
!test 05AUG2012    sgg= - (fnc_val - restmass_sph)
    sgg= - (restmass - restmass_sph)
    gg = (restmass - fnc_val)/hi
    ii = min0(5,count_adj)
    fac0 = facp(ii)
    fac0 = 1.0d0
!test 05AUG2012    emdc = par_val + sgg/gg*fac0
    emdc = emdc + sgg/gg*fac0
    epsmax_fnc = dmax1(epsmax_fnc,dabs(par_val - emdc)/dabs(par_val))
    if (epsmax_fnc.le.mass_eps) then 
      flag_restmass = 2
    else 
      flag_restmass = 0
    end if
!
    write(6,'(1a30,2(7x,i5))') ' -- Rest mass iteration # --  ',count_adj, ii
    write(6,'(1a30,1p,3e25.15)')' -- emdc   OLD -> NEW     --  ', par_val, emdc, dabs(par_val-emdc)
!
!write(6,*)flag_restmass, count_adj
!write(6,*)ii, fac0
!write(6,*)fixeddlm,fnc_val,hi
!write(6,*)restmass, restmass_sph
!write(6,*)admmass,  komarmass
!write(6,*)emdc,par_val
!write(6,*)sgg,gg
!write(6,*)epsmax_fnc, hi_fac
!!
    return
  end if
end subroutine adjust_rest_mass
