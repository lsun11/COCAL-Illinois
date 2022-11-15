subroutine printout_error_type2(iter_count,error_emd,char_type)
  implicit none
  real(8) :: error_emd
  integer :: iter_count
  character(len=3) :: char_type
!
!  write(6,'(a15,1p,e14.6)') ' Error       = ', error_emd
!  write(6,'(a15,i4)')       ' Iteration # = ', iter_count
  if (char_type.eq.'emd') then 
    write(6,'(a19,1p,e14.6)') ' Error in p/rho  = ', error_emd
  else
    write(6,'(a19,1p,e14.6)') ' Error           = ', error_emd
  end if
  write(6,'(a19,i4)')         ' Iteration #     = ', iter_count
end subroutine printout_error_type2
