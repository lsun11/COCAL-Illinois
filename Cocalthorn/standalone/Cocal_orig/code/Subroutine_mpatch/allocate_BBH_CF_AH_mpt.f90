subroutine allocate_BBH_CF_AH_mpt
  use phys_constant, only : nmpt
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_mpt
  use def_metric_excurve_grid_mpt
  use make_array_4d
  use make_array_6d
  implicit none
!
  call alloc_array4d(psi_ , 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(alph_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(alps_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvxd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvyd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvzd_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvxu_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvyu_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
  call alloc_array4d(bvzu_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!
!!!  call alloc_array4d(tfkijkij_, 1,nrg, 1,ntg, 1,npg, 1,nmpt)
!!!  call alloc_array6d(tfkij_,    1,nrg, 1,ntg, 1,npg, 1,3, 1,3, 1,nmpt)
!!!  call alloc_array4d(tfkijkij_grid_, 0,nrg, 0,ntg, 0,npg, 1,nmpt)
!!!  call alloc_array6d(tfkij_grid_, 0,nrg, 0,ntg, 0,npg, 1,3, 1,3, 1,nmpt)
!
end subroutine allocate_BBH_CF_AH_mpt
