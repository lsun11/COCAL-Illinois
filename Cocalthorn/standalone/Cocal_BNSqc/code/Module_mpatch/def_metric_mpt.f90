module def_metric_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  alph_(:,:,:,:), psi_(:,:,:,:), alps_(:,:,:,:)
  real(long), pointer  ::  alps2_(:,:,:,:)
  real(long), pointer  ::  bvxd_(:,:,:,:), bvyd_(:,:,:,:), bvzd_(:,:,:,:)
  real(long), pointer  ::  bvxu_(:,:,:,:), bvyu_(:,:,:,:), bvzu_(:,:,:,:)
  real(long), pointer  ::  tfkij_(:,:,:,:,:,:), tfkijkij_(:,:,:,:)
  real(long), pointer  ::  trk_(:,:,:,:)
end module def_metric_mpt
