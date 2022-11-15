! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_int_array_static_0dto1d_mpt
  use phys_constant, only : long, nnmpt
  implicit none
contains
! - - - - -
! 
subroutine copy_int_arraystatic_0dto1d_mpt(impt,array1,array2)
  implicit none
  integer, intent(in)  :: impt
  integer              :: array1, array2(1:nnmpt)
      array2(impt) &
  & = array1
end subroutine copy_int_arraystatic_0dto1d_mpt
end module copy_int_array_static_0dto1d_mpt
