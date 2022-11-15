subroutine printout_error(iter_count,error_emd)
  implicit none
  real(8) :: error_emd
  integer :: iter_count
!
!  write(6,'(a19,1p,e14.6)') ' Error       = ', error_emd
!  write(6,'(a19,i4)')       ' Iteration # = ', iter_count
  write(6,'(a20,i4,1p,e14.6)') ' Iteration, Error = ', iter_count, error_emd
end subroutine printout_error
