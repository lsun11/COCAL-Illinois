! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_int_array_static_1dto0d_mpt
  use phys_constant, only : long, nnmpt
  implicit none
contains
! - - - - -
! 
subroutine copy_int_arraystatic_1dto0d_mpt(impt,array1,array2)
  implicit none
  integer, intent(in)  :: impt
  integer              :: array1(1:nnmpt), array2
      array2 &
  & = array1(impt)
end subroutine copy_int_arraystatic_1dto0d_mpt
end module copy_int_array_static_1dto0d_mpt
