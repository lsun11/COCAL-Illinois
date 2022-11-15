! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_int_array_0dto1d_mpt
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 
subroutine copy_int_array0dto1d_mpt(impt,array1,array2)
  implicit none
  integer, intent(in)  :: impt
  integer, pointer  :: array1, array2(:)
      array2(impt) &
  & = array1
end subroutine copy_int_array0dto1d_mpt
end module copy_int_array_0dto1d_mpt
