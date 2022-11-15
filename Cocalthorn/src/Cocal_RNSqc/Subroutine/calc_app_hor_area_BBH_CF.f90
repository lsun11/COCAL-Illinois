subroutine calc_app_hor_area_BBH_CF
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg, rgin
  use make_array_2d
  use def_quantities, only : admmass, app_hor_area_bh, irredmass, bindingene
!  use coordinate_grav_r, only : rg
  use def_metric, only : psi
  use interface_interpo_linear_type0_2Dsurf
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: surf_int_bh, val
  real(long), pointer :: sou_bhsurf(:,:), psi_bh(:,:)
  integer :: irg, ipg, itg
!
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(psi_bh,0,ntg,0,npg)
!
  psi_bh(0:ntg,0:npg)  = psi(0,0:ntg,0:npg)
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,psi_bh,itg,ipg)
      sou_bhsurf(itg,ipg) = val**4
    end do
  end do
!
  call surf_int_grav_rg(sou_bhsurf,surf_int_bh,0)
  app_hor_area_bh = surf_int_bh
!
  irredmass  = 2.0d0*sqrt(app_hor_area_bh/16.0d0/pi)
  bindingene = admmass - irredmass
!
  deallocate(sou_bhsurf)
  deallocate(psi_bh)
end subroutine calc_app_hor_area_BBH_CF
