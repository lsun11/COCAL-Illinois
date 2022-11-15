subroutine allocate_pBH_CF
  use phys_constant, only : long
  use grid_parameter
  use def_metric_pBH
  use make_array_3d
  use make_array_5d
  implicit none
!
  call alloc_array3d(wme ,  0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(log_wme,  0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(log_N,  0, nrg, 0, ntg, 0, npg)
  call alloc_array5d(aij_trpBH,1,nrg,1,ntg,1,npg,1,3,1,3)
  call alloc_array3d(aijaij_trpBH,1,nrg,1,ntg,1,npg)
  call alloc_array5d(aij_trpBH_grid,0,nrg,0,ntg,0,npg,1,3,1,3)
  call alloc_array3d(aijaij_trpBH_grid,0,nrg,0,ntg,0,npg)
!
  call allocate_horizon
!
end subroutine allocate_pBH_CF
