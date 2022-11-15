subroutine allocate_SEM_tensor
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_SEM_tensor
  use make_array_3d
  use make_array_4d
  use make_array_5d
  implicit none
!
  call alloc_array3d(rhoH,  0, nrf, 0, ntf, 0, npf)
  call alloc_array4d(jmd,   0, nrf, 0, ntf, 0, npf, 1, 3)
  call alloc_array4d(jmu,   0, nrf, 0, ntf, 0, npf, 1, 3)
  call alloc_array5d(smijd, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, 3)
  call alloc_array5d(smiju, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, 3)
  call alloc_array3d(trsm,  0, nrf, 0, ntf, 0, npf)
!
!  call alloc_array3d(rhoH_grid,  0, nrf, 0, ntf, 0, npf)
!  call alloc_array3d(jmd_grid,   0, nrf, 0, ntf, 0, npf, 1, 3)
!  call alloc_array4d(jmu_grid,   0, nrf, 0, ntf, 0, npf, 1, 3)
!  call alloc_array4d(smijd_grid, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, 3)
!  call alloc_array5d(smiju_grid, 0, nrf, 0, ntf, 0, npf, 1, 3, 1, 3)
!  call alloc_array3d(trsm_grid,  0, nrf, 0, ntf, 0, npf)
!
end subroutine allocate_SEM_tensor
