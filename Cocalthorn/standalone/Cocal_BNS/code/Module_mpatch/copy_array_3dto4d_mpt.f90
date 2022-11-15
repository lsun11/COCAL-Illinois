! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_array_3dto4d_mpt
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 
subroutine copy_array3dto4d_mpt(impt,array1,array2,n1min,n1max,n2min,n2max,n3min,n3max)
  implicit none
  integer, intent(in)  :: n1min, n1max, n2min, n2max, n3min, n3max, impt
  real(long), pointer  :: array1(:,:,:), array2(:,:,:,:)
      array2(n1min:n1max,n2min:n2max,n3min:n3max,impt) &
  & = array1(n1min:n1max,n2min:n2max,n3min:n3max)
end subroutine copy_array3dto4d_mpt
end module copy_array_3dto4d_mpt
