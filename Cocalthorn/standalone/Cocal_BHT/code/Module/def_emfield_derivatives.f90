module def_emfield_derivatives
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  Lie_bAxd(:,:,:), Lie_bAyd(:,:,:), Lie_bAzd(:,:,:)
  real(long), pointer  ::  Lie_bAxd_grid(:,:,:), gLie_bAxu_grid(:,:,:), &
  &                        Lie_bAyd_grid(:,:,:), gLie_bAyu_grid(:,:,:), &
  &                        Lie_bAzd_grid(:,:,:), gLie_bAzu_grid(:,:,:)
  real(long), pointer  ::  cdvaxd(:,:,:,:), cdvayd(:,:,:,:), cdvazd(:,:,:,:)
  real(long), pointer  ::  pdvaxd(:,:,:,:), pdvayd(:,:,:,:), pdvazd(:,:,:,:)
end module def_emfield_derivatives
