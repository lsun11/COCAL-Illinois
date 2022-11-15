subroutine allocate_SEM_tensor_EMF
  use phys_constant, only : long
  use grid_parameter
  use def_SEM_tensor_EMF
  use make_array_3d
  use make_array_4d
  use make_array_5d
  implicit none
!
  call alloc_array3d(rhoH_EMF,  0, nrg, 0, ntg, 0, npg)
  call alloc_array4d(jmd_EMF,   0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array4d(jmu_EMF,   0, nrg, 0, ntg, 0, npg, 1, 3)
  call alloc_array5d(smijd_EMF, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 3)
  call alloc_array5d(smiju_EMF, 1, nrg, 1, ntg, 1, npg, 1, 3, 1, 3)
  call alloc_array3d(trsm_EMF,  0, nrg, 0, ntg, 0, npg)
!
!  call alloc_array3d(rhoH_EMF_grid,  0, nrg, 0, ntg, 0, npg)
!  call alloc_array3d(jmd_EMF_grid,   0, nrg, 0, ntg, 0, npg, 1, 3)
!  call alloc_array4d(jmu_EMF_grid,   0, nrg, 0, ntg, 0, npg, 1, 3)
!  call alloc_array4d(smijd_EMF_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
!  call alloc_array5d(smiju_EMF_grid, 0, nrg, 0, ntg, 0, npg, 1, 3, 1, 3)
!  call alloc_array3d(trsm_EMF_grid,  0, nrg, 0, ntg, 0, npg)
!
end subroutine allocate_SEM_tensor_EMF
