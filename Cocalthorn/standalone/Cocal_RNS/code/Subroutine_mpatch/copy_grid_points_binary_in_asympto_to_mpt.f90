subroutine copy_grid_points_binary_in_asympto_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_binary_in_asympto
  use grid_points_binary_in_asympto_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt,rb_a,  rb_a_,  -2,nrg+2,0,ntg,0,npg)
  call copy_array3dto4d_mpt(impt,thb_a, thb_a_, -2,nrg+2,0,ntg,0,npg)
  call copy_array3dto4d_mpt(impt,phib_a,phib_a_,-2,nrg+2,0,ntg,0,npg)
  call copy_array3dto4d_mpt(impt,hrb_a,  hrb_a_,  1,nrg,1,ntg,1,npg)
  call copy_array3dto4d_mpt(impt,hthb_a, hthb_a_, 1,nrg,1,ntg,1,npg)
  call copy_array3dto4d_mpt(impt,hphib_a,hphib_a_,1,nrg,1,ntg,1,npg)
!
end subroutine copy_grid_points_binary_in_asympto_to_mpt
