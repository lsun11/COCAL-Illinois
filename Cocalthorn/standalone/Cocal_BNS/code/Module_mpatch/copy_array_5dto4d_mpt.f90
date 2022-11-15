! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_array_5dto4d_mpt
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 
subroutine copy_array5dto4d_mpt(impt,array1,array2,n1min,n1max,n2min,n2max,n3min,n3max,n4min,n4max)
  implicit none
  integer, intent(in)  :: n1min, n1max, n2min, n2max, n3min, n3max, n4min, n4max, impt
  real(long), pointer  :: array1(:,:,:,:,:), array2(:,:,:,:)
      array2(n1min:n1max,n2min:n2max,n3min:n3max,n4min:n4max) &
  & = array1(n1min:n1max,n2min:n2max,n3min:n3max,n4min:n4max,impt)
end subroutine copy_array5dto4d_mpt
end module copy_array_5dto4d_mpt
