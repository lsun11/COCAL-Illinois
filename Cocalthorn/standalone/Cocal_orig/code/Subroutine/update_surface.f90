subroutine update_surface(potrs,rsnew,convf)
  use phys_constant, only : long
  use grid_parameter, only : ntg, npg
  implicit none
  real(long), pointer    :: potrs(:,:)
  real(long), pointer    :: rsnew(:,:)
  real(long), intent(in) :: convf
  rsnew(0:ntg,0:npg) =       convf *potrs(0:ntg,0:npg) &
  &                  + (1.d0-convf)*rsnew(0:ntg,0:npg)
end subroutine update_surface
