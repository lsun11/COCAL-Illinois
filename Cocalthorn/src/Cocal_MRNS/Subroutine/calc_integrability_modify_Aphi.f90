subroutine calc_integrability_modify_Aphi(vaxd_mod,vayd_mod,supp)
  use phys_constant, only  : long
  use grid_parameter
  use coordinate_grav_r, only : rg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_matter, only : rs
  use def_matter_parameter, only : ome
  use def_metric, only : alph, psi, bvxu, bvyu, bvzu
  use def_emfield, only : alva, vaxd, vayd, vazd
  use def_vector_phi, only : vec_phig
  use integrability_fnc_MHD
  use make_array_2d
  use make_array_3d
  implicit none
  real(long), pointer :: vaxd_mod(:,:,:), vayd_mod(:,:,:)
  real(long), pointer :: vaxd_0(:,:),     vayd_0(:,:)
  real(long) :: alvagg, bvxugg, bvyugg, bvzugg, vaxdgg, vazdgg, vphigg
  integer    :: irg, itg, ipg
  character(len=2), intent(in) :: supp
!
  if (MHDidx_q.ne.0.0d0) stop ' stop MHDidx_q '
!
  call alloc_array2d(vaxd_0, 0, nrg, 0, ntg)
  call alloc_array2d(vayd_0, 0, nrg, 0, ntg)
!
  vaxd_0(0:nrg,0:ntg) = vaxd(0:nrg,0:ntg,0)
  vayd_0(0:nrg,0:ntg) = vayd(0:nrg,0:ntg,0)
!
! compute va on phi=0 plane from integrability condition
! ### Can be used only for the case with A_t = -Omega Aphi + C
  ipg = 0
  do itg = 0, ntg
    do irg = 1, nrg
!
      if (supp.eq.'ns'.and.rg(irg).gt.rs(itg,ipg)) cycle
!
      alvagg = alva(irg,itg,ipg)
      bvxugg = bvxu(irg,itg,ipg)
      bvyugg = bvyu(irg,itg,ipg)
      bvzugg = bvzu(irg,itg,ipg)
      vaxdgg = vaxd(irg,itg,ipg)
      vazdgg = vazd(irg,itg,ipg)
      vphigg = vec_phig(irg,itg,ipg,2)
      vayd_0(irg,itg) = (alvagg - vaxdgg*bvxugg - vazdgg*bvzugg &
      &               +  MHDpar_charge)/(ome*vphigg + bvyugg)
!
    end do
  end do
!
! Copy to all phi planes.
!
  do ipg = 0, npg
!    if (supp.eq.'ns'.and.rg(irg).gt.rs(itg,ipg)) cycle
    vaxd_mod(0:nrg,0:ntg,ipg) = cosphig(ipg)*vaxd_0(0:nrg,0:ntg) &
    &                         - sinphig(ipg)*vayd_0(0:nrg,0:ntg)
    vayd_mod(0:nrg,0:ntg,ipg) = sinphig(ipg)*vaxd_0(0:nrg,0:ntg) &
    &                         + cosphig(ipg)*vayd_0(0:nrg,0:ntg)
  end do
!
  deallocate(vaxd_0)
  deallocate(vayd_0)
end subroutine calc_integrability_modify_Aphi
