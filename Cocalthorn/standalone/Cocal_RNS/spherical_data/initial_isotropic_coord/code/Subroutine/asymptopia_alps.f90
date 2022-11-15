subroutine asymptopia_alps(soupsou,souapou,dsoupsou,dsouapou)
!
  implicit none
  real(8), intent(out) :: soupsou, souapou, dsoupsou, dsouapou
!
! --- Compute source for GR boundary terms.
!
  soupsou = 1.0d0
  souapou = 1.0d0
  dsoupsou = 0.0d0
  dsouapou = 0.0d0
!
end subroutine asymptopia_alps
