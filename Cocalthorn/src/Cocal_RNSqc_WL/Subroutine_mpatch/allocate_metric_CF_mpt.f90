subroutine allocate_metric_CF_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_mpt
  use make_array_4d
  implicit none
!
  call alloc_array4d(psi_ , 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(alph_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(alps_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvxd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvyd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvzd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!
end subroutine allocate_metric_CF_mpt
