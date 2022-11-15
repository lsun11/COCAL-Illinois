subroutine allocate_BNS_CF_irrot
  use phys_constant, only : long
  use grid_parameter
  use def_matter
  use def_velocity_potential
!  use def_velocity_rot
  use make_array_2d
  use make_array_3d
  use make_array_4d
  implicit none
!
!  call alloc_array4d(vrot,0, nrf, 0, ntf, 0, npf, 1, 3)

  call alloc_array3d(vep  , 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepxf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepyf, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vepzf, 0, nrf, 0, ntf, 0, npf)

  call alloc_array3d(vepxg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vepyg, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(vepzg, 0, nrg, 0, ntg, 0, npg)

  call alloc_array2d(alm, 0, nlg, 0, nlg)
!
!  vrot(0:nrf,0:ntf,0:npf,1:3) = 0.0d0

end subroutine allocate_BNS_CF_irrot
