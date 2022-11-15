! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_int_array_1dto0d_mpt
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 
subroutine copy_int_array1dto0d_mpt(impt,array1,array2)
  implicit none
  integer, intent(in)  :: impt
  integer, pointer     :: array1(:), array2
      array2 &
  & = array1(impt)
end subroutine copy_int_array1dto0d_mpt
end module copy_int_array_1dto0d_mpt
