subroutine allocate_metric_3p1_WL_CTT
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric
  use def_metric_hij
  use def_metric_excurve_grid
  use def_CTT_decomposition, only : wvxd, wvyd, wvzd, sigt
  use make_array_3d
  use make_array_5d
  implicit none
!
  call alloc_array3d(psi, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hxxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hxyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hxzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hyyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hyzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hzzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hxxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hxyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hxzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hyyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hyzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(hzzu, 0, nrg, 0, ntg, 0, npg)
!
  call alloc_array3d(gaugex, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gaugey, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gaugez, 0, nrg, 0, ntg, 0, npg)

  call alloc_array3d(wvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(wvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(wvzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(sigt, 0, nrg, 0, ntg, 0, npg)

  call alloc_array5d(tfkij_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
  call alloc_array3d(trk_grid, 0, nrg, 0, ntg, 0, npg)
!
end subroutine allocate_metric_3p1_WL_CTT
