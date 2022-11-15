subroutine torus_emd_max_c
  use phys_constant,  only  : long
  use grid_parameter, only  : ntgxy,nrg,nrgin
  use def_matter, only  :   emdg, omeg, jomeg
  use def_matter_parameter, only : emdc
  use def_bht_parameter
  use coordinate_grav_r, only : rg  
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  implicit none
  integer    :: irem, ithem, iphem, nmx, irg
!
  o_in = omeg(nrgin,ntgxy,0)   
  j_in = jomeg(nrgin,ntgxy,0)
  l_in = j_in/(1.0d0+o_in*j_in)

  do irg = 0, nrg-2
    if (emdg(irg,ntgxy,0).gt.emdg(irg+1,ntgxy,0) .and. emdg(irg+1,ntgxy,0)==emdg(irg+2,ntgxy,0)) then
      o_out = omeg(irg+1,ntgxy,0)
      j_out = jomeg(irg+1,ntgxy,0)
      l_out = j_out/(1.0d0+o_out*j_out)
      exit
    end if
  end do

  call search_emdmax_xaxis_grgrid(nmx)
  xe_c = rg(nmx);  ye_c=0.0d0;   ze_c=0.0d0
  o_c = omeg(nmx,ntgxy,0)   
  j_c = jomeg(nmx,ntgxy,0)
  l_c = j_c/(1.0d0+o_c*j_c)

  call search_emdmax_full_grgrid(irem,ithem,iphem)
  xe_m = rg(irem)*sinthg(ithem)*cosphig(iphem)
  ye_m = rg(irem)*sinthg(ithem)*sinphig(iphem)
  ze_m = rg(irem)*costhg(ithem)
  o_m = omeg(irem,ithem,iphem) 
  j_m = jomeg(irem,ithem,iphem)
  l_m = j_m/(1.0d0+o_m*j_m)

!    emdc = emdg(nmx,ntgxy,0)   ! slightly different from emdc of calculated solution
  write(6,'(a24,1p3e23.15)') " Max emdg on x-axis   = ", emdg(nmx,ntgxy,0)
  write(6,'(a24,1p3e23.15)') " ...at xe_c,ye_c,ze_c = ", xe_c,ye_c,ze_c
  write(6,'(a24,1p3e23.15)') " emdc                 = ", emdc
  write(6,'(a24,1p3e23.15)') " Max emdg             = ", emdg(irem,ithem,iphem)
  write(6,'(a24,1p3e23.15)') " ...at xe_m,ye_m,ze_m = ", xe_m, ye_m, ze_m
!
end subroutine torus_emd_max_c
