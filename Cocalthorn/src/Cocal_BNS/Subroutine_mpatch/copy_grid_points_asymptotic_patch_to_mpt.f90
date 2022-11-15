subroutine copy_grid_points_asymptotic_patch_to_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_asymptotic_patch
  use grid_points_asymptotic_patch_mpt
  use copy_array_3dto4d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt,ra,  ra_,  -2,nrg+2,0,ntg,0,npg)
  call copy_array3dto4d_mpt(impt,tha, tha_, -2,nrg+2,0,ntg,0,npg)
  call copy_array3dto4d_mpt(impt,phia,phia_,-2,nrg+2,0,ntg,0,npg)
  call copy_array3dto4d_mpt(impt,hra,  hra_,  1,nrg,1,ntg,1,npg)
  call copy_array3dto4d_mpt(impt,htha, htha_, 1,nrg,1,ntg,1,npg)
  call copy_array3dto4d_mpt(impt,hphia,hphia_,1,nrg,1,ntg,1,npg)
!
end subroutine copy_grid_points_asymptotic_patch_to_mpt
