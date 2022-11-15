subroutine printout_error_metric_combined_3(iter_count,ebx,eby,ebz)
  implicit none
  real(8) :: ebx,eby,ebz
  integer :: iter_count
  write(6,'(a14,i5,1p,3e15.6)') ' Iteration # =', iter_count, ebx, eby, ebz
end subroutine printout_error_metric_combined_3
