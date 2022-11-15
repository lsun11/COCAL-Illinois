!  Associated Legendre function and factorials
!______________________________________________
subroutine copy_legendre_fn_grav_from_mpt(impt)
  use grid_parameter, only : nlg, ntg
  use legendre_fn_grav
  use legendre_fn_grav_mpt
  use copy_array_4dto3d_mpt
  use copy_array_2dto1d_mpt
  use copy_array_3dto2d_mpt
  implicit none
  integer :: impt
!
  call copy_array4dto3d_mpt(impt, plmg_, plmg, 0, nlg, 0, nlg, 0, ntg)
  call copy_array4dto3d_mpt(impt, hplmg_, hplmg, 0, nlg, 0, nlg, 1, ntg)
  call copy_array3dto2d_mpt(impt, facnmg_, facnmg, 0, nlg, 0, nlg)
  call copy_array2dto1d_mpt(impt, epsig_, epsig, 0, nlg)
!
end subroutine copy_legendre_fn_grav_from_mpt
