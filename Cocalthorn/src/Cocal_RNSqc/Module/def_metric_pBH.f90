module def_metric_pBH
  use phys_constant, only : long, nnrg
  implicit none
  real(long) ::    R_trpBH(0:nnrg),    hR_trpBH(nnrg)
  real(long) ::  psi_trpBH(0:nnrg),  hpsi_trpBH(nnrg)
  real(long) :: alph_trpBH(0:nnrg), halph_trpBH(nnrg)
  real(long) ::  bvr_trpBH(0:nnrg),  hbvr_trpBH(nnrg)
  real(long), pointer  ::  wme(:,:,:), log_wme(:,:,:), log_N(:,:,:)
  real(long), pointer  ::  aij_trpBH(:,:,:,:,:), aijaij_trpBH(:,:,:)
  real(long), pointer  ::  aij_trpBH_grid(:,:,:,:,:), aijaij_trpBH_grid(:,:,:)
  integer, parameter :: index_wme = 2
end module def_metric_pBH
