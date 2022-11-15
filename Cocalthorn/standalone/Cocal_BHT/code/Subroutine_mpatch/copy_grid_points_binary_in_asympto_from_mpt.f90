subroutine copy_grid_points_binary_in_asympto_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_binary_in_asympto
  use grid_points_binary_in_asympto_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt,rb_a_,  rb_a,  -2,nrg+2,0,ntg,0,npg)
  call copy_array4dto3d_mpt(impt,thb_a_, thb_a, -2,nrg+2,0,ntg,0,npg)
  call copy_array4dto3d_mpt(impt,phib_a_,phib_a,-2,nrg+2,0,ntg,0,npg)
  call copy_array4dto3d_mpt(impt,hrb_a_,  hrb_a,  1,nrg,1,ntg,1,npg)
  call copy_array4dto3d_mpt(impt,hthb_a_, hthb_a, 1,nrg,1,ntg,1,npg)
  call copy_array4dto3d_mpt(impt,hphib_a_,hphib_a,1,nrg,1,ntg,1,npg)
!
end subroutine copy_grid_points_binary_in_asympto_from_mpt
