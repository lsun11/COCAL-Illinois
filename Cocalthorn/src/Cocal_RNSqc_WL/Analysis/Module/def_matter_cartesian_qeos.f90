module def_matter_cartesian_qeos
  use phys_constant, only : long
  implicit none
  real(long), pointer  :: rhoca(:,:,:), utca(:,:,:), vepca(:,:,:)
  real(long), pointer  ::  vxca(:,:,:), vyca(:,:,:), vzca(:,:,:)
  real(long), pointer  ::  wxspca(:,:,:), wyspca(:,:,:), wzspca(:,:,:)
end module def_matter_cartesian_qeos
