! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_int_array_3dto2d_mpt
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 
subroutine copy_int_array3dto2d_mpt(impt,array1,array2,n1min,n1max,n2min,n2max)
  implicit none
  integer, intent(in)  :: n1min, n1max, n2min, n2max, impt
  integer, pointer  :: array1(:,:,:), array2(:,:)
      array2(n1min:n1max,n2min:n2max) &
  & = array1(n1min:n1max,n2min:n2max,impt)
end subroutine copy_int_array3dto2d_mpt
end module copy_int_array_3dto2d_mpt
