subroutine reset_outer_boundary_BBH_CF(char_mp)
  use phys_constant, only  : long, nmpt
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd
  use def_matter_parameter, only : radi
  use def_quantities, only : admmass_asymp, komarmass_asymp
  implicit none
  real(long), pointer :: sou_surf(:,:)
  character(len=4), intent(in) :: char_mp
! Reset boundary condition for BH
! ###NOT FOT NS!!!  For NS, rg should be replaced by radi*rg
!
  call copy_def_quantities_from_mpt(nmpt)
!
  if (char_mp.eq.'psi ') then
    psi(nrg,0:ntg,0:npg) = 1.0d0 + admmass_asymp/(2.0d0*rg(nrg))
  end if
  if (char_mp.eq.'alph') then
    alph(nrg,0:ntg,0:npg)=(1.0d0 - (2.0d0*komarmass_asymp-admmass_asymp) &
    &                             /(2.0d0*rg(nrg)))/psi(nrg,0:ntg,0:npg)
  end if
  if (char_mp.eq.'bvxd') bvxd(nrg,0:ntg,0:npg) = 0.0d0
  if (char_mp.eq.'bvyd') bvyd(nrg,0:ntg,0:npg) = 0.0d0
  if (char_mp.eq.'bvzd') bvzd(nrg,0:ntg,0:npg) = 0.0d0
!
end subroutine reset_outer_boundary_BBH_CF
