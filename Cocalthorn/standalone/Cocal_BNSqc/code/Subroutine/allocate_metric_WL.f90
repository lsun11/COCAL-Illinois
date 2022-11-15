subroutine allocate_metric_WL
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric
  use def_gamma_crist
  use def_gamma_crist_grid
  use def_metric_hij
  use def_metric_hij_dirac
  use def_metric_rotshift
  use def_ricci_tensor
  use def_metric_excurve_grid
  use def_shift_derivatives
  use def_shift_derivatives_grid
  use def_Lie_derivatives
  use def_Lie_derivatives_grid
  use def_cristoffel
  use def_cristoffel_grid
  use def_transverse_part
  use make_array_2d
  use make_array_3d
  use make_array_4d
  use make_array_5d
  implicit none
!
  call alloc_array3d(psi, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alph, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(alps2, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(bvzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(ovxu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(ovyu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(ovzu, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(ovxd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(ovyd, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(ovzd, 0, nrg, 0, ntg, 0, npg)
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
!
  call alloc_array3d(tfkijkij, 1, nrg, 1, ntg, 1, npg)
  call alloc_array5d(tfkij, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 3)
  call alloc_array3d(tfkijkij_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array5d(tfkij_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
  call alloc_array3d(trk, 1, nrg, 1, ntg, 1, npg)
  call alloc_array3d(trk_grid, 0, nrg, 0, ntg, 0, npg)
!
  call alloc_array5d(cri, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 6)
  call alloc_array5d(cri_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 6)
  call alloc_array5d(crid, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 6)
  call alloc_array5d(crid_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 6)
  call alloc_array3d(gmcrix, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gmcriy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gmcriz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gmcrix_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gmcriy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(gmcriz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array4d(dagmabu, 1, nrg, 1, ntg, 1, npg, 1, 3)
!
  call alloc_array3d(Ftvx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Ftvy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Ftvz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Ftvx_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Ftvy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(Ftvz_grid, 0, nrg, 0, ntg, 0, npg)
!
  call alloc_array4d(rab, 0, nrg, 0, ntg, 0, npg, 1, 6)
  call alloc_array4d(rabnl, 0, nrg, 0, ntg, 0, npg, 1, 6)
  call alloc_array4d(rabDF, 0, nrg, 0, ntg, 0, npg, 1, 6)
  call alloc_array3d(elpxx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpxy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpxz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpyy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpyz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpzz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpxx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpxy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpxz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpyy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpyz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpzz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbxx, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbxy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbxz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbyy, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbyz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbzz, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpxx_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpxy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpxz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpyy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpyz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(elpzz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpxx_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpxy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpxz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpyy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpyz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlpzz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbxx_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbxy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbxz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbyy_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbyz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(rlbzz_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array4d(cdbvxd, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(cdbvyd, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(cdbvzd, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(cdbvxd_grid, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(cdbvyd_grid, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(cdbvzd_grid, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array3d(cdivbv, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(cdivbv_grid, 0, nrg, 0, ntg, 0, npg)
  call alloc_array4d(pdbvxd_grid, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(pdbvyd_grid, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(pdbvzd_grid, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(pdbvxd, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(pdbvyd, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call alloc_array4d(pdbvzd, 1, nrg, 1, ntg, 1, npg, 1, 3)
!
end subroutine allocate_metric_WL
