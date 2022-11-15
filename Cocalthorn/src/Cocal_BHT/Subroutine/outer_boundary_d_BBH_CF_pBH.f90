subroutine outer_boundary_d_BBH_CF_pBH(sou_surf,char_mp)
  use phys_constant, only  : long, nmpt
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_matter_parameter, only : radi
  use def_metric_pBH, only : index_wme
  use def_quantities, only : admmass_asymp, komarmass_asymp
  implicit none
  real(long), pointer :: sou_surf(:,:)
  real(long) :: index, fac2
  character(len=4), intent(in) :: char_mp
! Reset boundary condition for BH
! ###NOT FOT NS!!!  For NS, rg should be replaced by radi*rg
!
  call copy_def_quantities_from_mpt(nmpt)
!
  index = dble(index_wme)
  fac2  = 1.0d0/dsqrt(2.0d0)
  if (char_mp.eq.'logw') then
    sou_surf(1:ntg,1:npg) = - index*admmass_asymp/(2.0d0*rg(nrg))
  end if
  if (char_mp.eq.'logN') then
    sou_surf(1:ntg,1:npg) = - fac2*komarmass_asymp/rg(nrg)
  end if
!
end subroutine outer_boundary_d_BBH_CF_pBH
