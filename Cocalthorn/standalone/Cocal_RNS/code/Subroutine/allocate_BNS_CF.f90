subroutine allocate_BNS_CF
  use phys_constant, only : long
  use grid_parameter
  use def_metric
  use def_matter
  use def_metric_excurve_grid
  use def_vector_x
  use def_vector_phi
!  use def_vector_bh
  use def_vector_irg
  use make_array_2d
  use make_array_3d
  use make_array_5d
  implicit none
!
  call alloc_array2d(rs,        0, ntf, 0, npf)
  call alloc_array3d(emdg,      0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(emd,       0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(utg,       0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(utf,       0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(omeg,      0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(omef,      0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg,     0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef,     0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(jomeg_int, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(jomef_int, 0, nrf, 0, ntf, 0, npf)
!
  call alloc_array3d(psi , 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzu, 0, nrg, 0, ntg, 0, npg)
!
  call alloc_array3d(tfkijkij,1,nrg,1,ntg,1,npg)
  call alloc_array5d(tfkij,   1,nrg,1,ntg,1,npg,1,3,1,3)
  call alloc_array3d(tfkijkij_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array5d(tfkij_grid,    0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
  call alloc_array3d(trk,1,nrg,1,ntg,1,npg)
  call alloc_array3d(trk_grid, 0, nrg, 0, ntg, 0, npg)
!  
  call allocate_vector_x
  call allocate_vector_phi
!  call allocate_vector_bh
  call allocate_vector_irg
!
  psi(0:nrg,0:ntg,0:npg) =1.0d0
  alph(0:nrg,0:ntg,0:npg)=1.0d0
  alps(0:nrg,0:ntg,0:npg)=1.0d0
  bvxd(0:nrg,0:ntg,0:npg)=0.0d0
  bvyd(0:nrg,0:ntg,0:npg)=0.0d0
  bvzd(0:nrg,0:ntg,0:npg)=0.0d0

  trk(1:nrg,1:ntg,1:npg)=0.0d0
  tfkijkij(1:nrg,1:ntg,1:npg)     =0.0d0
  tfkij(1:nrg,1:ntg,1:npg,1:3,1:3)=0.0d0

  trk_grid(0:nrg,0:ntg,0:npg)=0.0d0
  tfkijkij_grid(0:nrg,0:ntg,0:npg)     =0.0d0
  tfkij_grid(0:nrg,0:ntg,0:npg,1:3,1:3)=0.0d0

end subroutine allocate_BNS_CF
