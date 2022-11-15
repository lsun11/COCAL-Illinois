subroutine printout_error_metric_combined(iter_count,ep,ea,ebx,eby,ebz)
  implicit none
  real(8) :: ep,ea,ebx,eby,ebz
  integer :: iter_count
!  write(6,'(a12,i4,a10,a17,1p,e14.6)') 'Iteration #=', iter_count, '         ', 'Error in metric= ', error_pot
  write(6,'(a12,i3,1p,5e15.6)') 'Iteration #=', iter_count, ep,ea,ebx,eby,ebz
end subroutine printout_error_metric_combined
