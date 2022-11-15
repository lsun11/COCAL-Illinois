!  Associated Legendre function and factorials
!______________________________________________
subroutine copy_legendre_fn_grav_to_mpt(impt)
  use grid_parameter, only : nlg, ntg
  use legendre_fn_grav
  use legendre_fn_grav_mpt
  use copy_array_3dto4d_mpt
  use copy_array_1dto2d_mpt
  use copy_array_2dto3d_mpt
  implicit none
  integer :: impt
!
  call copy_array3dto4d_mpt(impt, plmg, plmg_, 0, nlg, 0, nlg, 0, ntg)
  call copy_array3dto4d_mpt(impt, hplmg, hplmg_, 0, nlg, 0, nlg, 1, ntg)
  call copy_array2dto3d_mpt(impt, facnmg, facnmg_, 0, nlg, 0, nlg)
  call copy_array1dto2d_mpt(impt, epsig, epsig_, 0, nlg)
!
end subroutine copy_legendre_fn_grav_to_mpt
