module grid_temporary
  use phys_constant, only : long, nnrg, nntg, nnpg
  implicit none
  integer :: nrgtmp, ntgtmp, npgtmp, nrftmp, ntftmp, npftmp
  real(long) :: rgtmp(0:nnrg), thgtmp(0:nntg), phigtmp(0:nnpg)
end module grid_temporary
