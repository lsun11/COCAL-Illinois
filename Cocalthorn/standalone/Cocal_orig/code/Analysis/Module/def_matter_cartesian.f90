module def_matter_cartesian
  use phys_constant, only : long
  implicit none
  real(long), pointer  :: emdca(:,:,:), utca(:,:,:), vepca(:,:,:)
  real(long), pointer  ::  vxca(:,:,:), vyca(:,:,:), vzca(:,:,:)
  real(long), pointer  ::  wxspca(:,:,:), wyspca(:,:,:), wzspca(:,:,:)
  real(long), pointer  :: omeca(:,:,:)
end module def_matter_cartesian
