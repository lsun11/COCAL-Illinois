subroutine calc_integrability_modify_At(va_mod,supp)
  use phys_constant, only  : long
  use grid_parameter
  use coordinate_grav_r, only : rg
  use def_matter, only : rs
  use def_metric, only : alph, psi, bvxu, bvyu, bvzu
  use def_emfield, only : va, vaxd, vayd, vazd
  use def_vector_phi, only : vec_phig
  use integrability_fnc_MHD
  use make_array_3d
  implicit none
  real(long), pointer :: va_mod(:,:,:)
  real(long), pointer :: vaphid(:,:,:)
  real(long) :: Aphi
  real(long) :: alphgg, bvxugg, bvyugg, bvzugg, vaxdgg, vaydgg, vazdgg
  integer    :: irg, itg, ipg
  integer    :: irg1
  character(len=2), intent(in) :: supp
!
  call alloc_array3d(vaphid, 0, nrg, 0, ntg, 0, npg)
!
  vaphid(0:nrg,0:ntg,0:npg)=vayd(0:nrg,0:ntg,0:npg) &
  &                    *vec_phig(0:nrg,0:ntg,0:npg,2)
!
!do not use  va_mod(0:nrg,0:ntg,0:npg) = va(0:nrg,0:ntg,0:npg)
!
! compute va on phi=0 plane from integrability condition
  ipg = 0
  do itg = 0, ntg
    do irg = 0, nrg
!
      if (supp.eq.'ns'.and.rg(irg).gt.rs(itg,ipg)) cycle
!
      alphgg = alph(irg,itg,ipg)
      bvxugg = bvxu(irg,itg,ipg)
      bvyugg = bvyu(irg,itg,ipg)
      bvzugg = bvzu(irg,itg,ipg)
      vaxdgg = vaxd(irg,itg,ipg)
      vaydgg = vayd(irg,itg,ipg)
      vazdgg = vazd(irg,itg,ipg)
      Aphi = vaphid(irg,itg,ipg)
      call calc_integrability_fnc_MHD(Aphi)
      va_mod(irg,itg,ipg) = 1.0/alphgg*(- MHDfnc_At &
      &                   + vaxdgg*bvxugg + vaydgg*bvyugg + vazdgg*bvzugg)
!
    end do
  end do
!
! Copy to phi /= 0 planes.
!
  do ipg = 1, npg
    if (supp.eq.'ns'.and.rg(irg).gt.rs(itg,ipg)) cycle
    va_mod(0:nrg,0:ntg,ipg) = va_mod(0:nrg,0:ntg,0)
  end do
!
  deallocate(vaphid)
!
end subroutine calc_integrability_modify_At
