subroutine copy_def_vector_x_from_mpt(impt)
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_vector_x
  use def_vector_x_mpt
  use copy_array_5dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array5dto4d_mpt(impt, vec_xg_, vec_xg, 0, nrg, 0, ntg, 0, npg, 1, 3)
  call copy_array5dto4d_mpt(impt, vec_xf_, vec_xf, 0, nrf, 0, ntf, 0, npf, 1, 3)
  call copy_array5dto4d_mpt(impt, hvec_xg_, hvec_xg, 1, nrg, 1, ntg, 1, npg, 1, 3)
  call copy_array5dto4d_mpt(impt, hvec_xf_, hvec_xf, 1, nrf, 1, ntf, 1, npf, 1, 3)
!
end subroutine copy_def_vector_x_from_mpt
