subroutine allocate_velocity_potential
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nlg
  use def_velocity_potential
  use make_array_2d
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(vep, 0,nrf, 0,ntf, 0,npf)
  call alloc_array4d(grad_vep, 0,nrf, 0,ntf, 0,npf, 1,3)
  call alloc_array2d(alm, 0,nlg, 0,nlg)
!
end subroutine allocate_velocity_potential
