subroutine copy_grid_points_asymptotic_patch_from_mpt(impt)
  use grid_parameter, only : nrg, ntg, npg
  use grid_points_asymptotic_patch
  use grid_points_asymptotic_patch_mpt
  use copy_array_4dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt,ra_,  ra,  -2,nrg+2,0,ntg,0,npg)
  call copy_array4dto3d_mpt(impt,tha_, tha, -2,nrg+2,0,ntg,0,npg)
  call copy_array4dto3d_mpt(impt,phia_,phia,-2,nrg+2,0,ntg,0,npg)
  call copy_array4dto3d_mpt(impt,hra_,  hra,  1,nrg,1,ntg,1,npg)
  call copy_array4dto3d_mpt(impt,htha_, htha, 1,nrg,1,ntg,1,npg)
  call copy_array4dto3d_mpt(impt,hphia_,hphia,1,nrg,1,ntg,1,npg)
!
end subroutine copy_grid_points_asymptotic_patch_from_mpt
