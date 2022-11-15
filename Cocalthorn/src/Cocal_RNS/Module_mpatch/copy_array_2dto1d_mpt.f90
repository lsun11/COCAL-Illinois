! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_array_2dto1d_mpt
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 
subroutine copy_array2dto1d_mpt(impt,array1,array2,n1min,n1max)
  implicit none
  integer, intent(in)  :: n1min, n1max, impt
  real(long), pointer  :: array1(:,:), array2(:)
      array2(n1min:n1max) &
  & = array1(n1min:n1max,impt)
end subroutine copy_array2dto1d_mpt
end module copy_array_2dto1d_mpt
