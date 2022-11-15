! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_array_3d
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 3D array
subroutine copy_array3d(array1,array2,n1min,n1max,n2min,n2max,n3min,n3max)
  implicit none
  integer, intent(in)  :: n1min, n1max, n2min, n2max, n3min, n3max
  real(long), pointer  :: array1(:,:,:), array2(:,:,:)
      array2(n1min:n1max,n2min:n2max,n3min:n3max) &
  & = array1(n1min:n1max,n2min:n2max,n3min:n3max)
end subroutine copy_array3d
end module copy_array_3d
