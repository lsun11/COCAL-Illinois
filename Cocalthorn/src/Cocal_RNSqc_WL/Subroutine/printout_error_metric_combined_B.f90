subroutine printout_error_metric_combined_B(iter_count,ep,ea,ebx,eby,ebz,eB)
  implicit none
  real(8) :: ep,ea,ebx,eby,ebz,eB
  integer :: iter_count
!  write(6,'(a12,i4,a10,a17,1p,e14.6)') 'Iteration #=', iter_count, '         ', 'Error in metric= ', error_pot
  write(6,'(a12,i3,1p,6e15.6)') 'Iteration #=', iter_count, ep,ea,ebx,eby,ebz,eB
end subroutine printout_error_metric_combined_B
