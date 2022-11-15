! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_array_4d
  use phys_constant, only : long
  implicit none
contains
! - - - - -
! 3D array
subroutine copy_array4d(array1,array2,n1min,n1max,n2min,n2max,&
&                                     n3min,n3max,n4min,n4max)
  implicit none
  integer, intent(in)  :: n1min, n1max, n2min, n2max, &
  &                       n3min, n3max, n4min, n4max
  real(long), pointer  :: array1(:,:,:,:), array2(:,:,:,:)
      array2(n1min:n1max,n2min:n2max,n3min:n3max,n4min:n4max) &
  & = array1(n1min:n1max,n2min:n2max,n3min:n3max,n4min:n4max)
end subroutine copy_array4d
end module copy_array_4d
