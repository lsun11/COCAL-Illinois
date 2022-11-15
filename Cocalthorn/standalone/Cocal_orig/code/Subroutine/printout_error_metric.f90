subroutine printout_error_metric(iter_count,error_pot)
  implicit none
  real(8) :: error_pot
  integer :: iter_count
!  write(6,'(a19,i4)')       ' Iteration #     = ', iter_count
  write(6,'(a19,1p,e14.6)') ' Error in metric = ', error_pot
end subroutine printout_error_metric
