subroutine initialize_field
  implicit none
!
! --- Preparation for other index
!
  call invhij
  call calc_shift_down2up
  call calc_shift2rotshift
!
! --- Preparation for the scalar part
!
  call ap2alps
!
end subroutine initialize_field
