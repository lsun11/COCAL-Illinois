module def_metric_pBH_mpt
  use phys_constant, only : long, nnrg, nnmpt
  implicit none
  real(long) ::    R_trpBH_(0:nnrg,nnmpt),    hR_trpBH_(nnrg,nnmpt)
  real(long) ::  psi_trpBH_(0:nnrg,nnmpt),  hpsi_trpBH_(nnrg,nnmpt)
  real(long) :: alph_trpBH_(0:nnrg,nnmpt), halph_trpBH_(nnrg,nnmpt)
  real(long) ::  bvr_trpBH_(0:nnrg,nnmpt),  hbvr_trpBH_(nnrg,nnmpt)
  real(long), pointer :: wme_(:,:,:,:), log_wme_(:,:,:,:), log_N_(:,:,:,:)
  real(long), pointer :: aij_trpBH_(:,:,:,:,:,:), aijaij_trpBH_(:,:,:,:)
  real(long), pointer :: aij_trpBH_grid_(:,:,:,:,:,:)
  real(long), pointer :: aijaij_trpBH_grid_(:,:,:,:)
end module def_metric_pBH_mpt
