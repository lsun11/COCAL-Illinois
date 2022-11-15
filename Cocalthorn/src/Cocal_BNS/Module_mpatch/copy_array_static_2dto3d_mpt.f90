! ----------------------------------------------------------------------
!  Copy array
! ----------------------------------------------------------------------
module copy_array_static_2dto3d_mpt
  use phys_constant, only : long, nnmpt
  implicit none
contains
! - - - - -
! 
subroutine copy_arraystatic_2dto3d_mpt(impt,array1,array2,n1min,n1max,n2min,n2max)
  implicit none
  integer, intent(in)  :: n1min, n1max, n2min, n2max, impt
  real(long) :: array1(n1min:n1max,n2min:n2max), &
  &             array2(n1min:n1max,n2min:n2max,1:nnmpt)
      array2(n1min:n1max,n2min:n2max,impt) &
  & = array1(n1min:n1max,n2min:n2max)
end subroutine copy_arraystatic_2dto3d_mpt
end module copy_array_static_2dto3d_mpt
