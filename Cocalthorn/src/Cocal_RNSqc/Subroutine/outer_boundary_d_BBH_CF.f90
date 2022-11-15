subroutine outer_boundary_d_BBH_CF(sou_surf,char_mp)
  use phys_constant, only  : long, nmpt
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
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
    sou_surf(1:ntg,1:npg) = 1.0d0 + admmass_asymp/(2.0d0*rg(nrg))
  end if
  if (char_mp.eq.'alps') then
    sou_surf(1:ntg,1:npg) = 1.0d0 - (2.0d0*komarmass_asymp-admmass_asymp) &
    &                                            /(2.0d0*rg(nrg))
  end if
  if (char_mp.eq.'bvxd'.or.char_mp.eq.'bvyd'.or.char_mp.eq.'bvzd') then
    sou_surf(1:ntg,1:npg) = 0.0d0
  end if
!
end subroutine outer_boundary_d_BBH_CF
