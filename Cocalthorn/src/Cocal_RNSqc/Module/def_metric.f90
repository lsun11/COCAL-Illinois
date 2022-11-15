module def_metric
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  alph(:,:,:), psi(:,:,:), alps(:,:,:)
  real(long), pointer  ::  alps2(:,:,:)
  real(long), pointer  ::  bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:)
  real(long), pointer  ::  bvxu(:,:,:), bvyu(:,:,:), bvzu(:,:,:)
  real(long), pointer  ::  tfkij(:,:,:,:,:), tfkijkij(:,:,:)
  real(long), pointer  ::  trk(:,:,:), ditrk(:,:,:,:)
end module def_metric
