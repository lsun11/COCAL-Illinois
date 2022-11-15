subroutine copy_def_vector_x_to_mpt(impt)
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_vector_x
  use def_vector_x_mpt
  use copy_array_4dto5d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto5d_mpt(impt, vec_xg, vec_xg_, 0,nrg,0,ntg,0,npg,1,3)
  call copy_array4dto5d_mpt(impt, vec_xf, vec_xf_, 0,nrf,0,ntf,0,npf,1,3)
  call copy_array4dto5d_mpt(impt,hvec_xg,hvec_xg_, 1,nrg,1,ntg,1,npg,1,3)
  call copy_array4dto5d_mpt(impt,hvec_xf,hvec_xf_, 1,nrf,1,ntf,1,npf,1,3)
!
end subroutine copy_def_vector_x_to_mpt
