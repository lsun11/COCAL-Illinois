module def_CTT_decomposition
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  sigt(:,:,:)
  real(long), pointer  ::  wvxd(:,:,:), wvyd(:,:,:), wvzd(:,:,:)
  real(long), pointer  ::  wvxu(:,:,:), wvyu(:,:,:), wvzu(:,:,:)
  real(long), pointer  ::  tfkij_CTT(:,:,:,:,:)
  real(long), pointer  ::  tfkij_grid_CTT(:,:,:,:,:)
end module def_CTT_decomposition
