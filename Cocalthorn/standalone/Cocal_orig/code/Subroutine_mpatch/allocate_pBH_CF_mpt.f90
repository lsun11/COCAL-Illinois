subroutine allocate_pBH_CF_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH_mpt
  use make_array_4d
  use make_array_6d
  implicit none
!
  call alloc_array4d(wme_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(log_wme_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(log_N_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!
!!!  call alloc_array4d(aijaij_trpBH_, 1,nrg, 1,ntg, 1,npg, 1,nmpt)
!!!  call alloc_array6d(aij_trpBH_,    1,nrg, 1,ntg, 1,npg, 1,3, 1,3, 1,nmpt)
!!!  call alloc_array4d(aijaij_trpBH_grid_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!!!  call alloc_array6d(aij_trpBH_grid_, 0,nrg, 0,ntg, 0,npg, 1,3, 1,3, 1,nmpt)
!
  call allocate_horizon_mpt
!
end subroutine allocate_pBH_CF_mpt
