subroutine outer_boundary_WL(sou_surf,dsou_surf,char_mp)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use def_binary_parameter, only : dis
  use def_metric
  use coordinate_grav_r, only : rg
  use def_quantities, only : admmass_asymp, komarmass_asymp
  implicit none
  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long), pointer :: psi_bhsurf(:,:), dpsi_bhsurf(:,:)
  character(len=4), intent(in) :: char_mp
!
  if (char_mp.eq.'psi ') then
    sou_surf(1:ntg,1:npg) = 1.0d0 + admmass_asymp/(2.0d0*rg(nrg))
!    sou_surf(1:ntg,1:npg) = 1.0d0
  else if (char_mp.eq.'alps') then
    sou_surf(1:ntg,1:npg) = 1.0d0 - (2.0d0*komarmass_asymp-admmass_asymp)/(2.0d0*rg(nrg))
  else   
    sou_surf(1:ntg,1:npg) = 0.0d0
  end if
!
end subroutine outer_boundary_WL
