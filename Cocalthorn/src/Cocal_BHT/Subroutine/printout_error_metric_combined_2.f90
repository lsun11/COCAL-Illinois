subroutine printout_error_metric_combined_2(iter_count,ep,ea)
  implicit none
  real(8) :: ep,ea
  integer :: iter_count
  write(6,'(a14,i5,1p,2e15.6)') ' Iteration # =', iter_count, ep, ea
end subroutine printout_error_metric_combined_2
