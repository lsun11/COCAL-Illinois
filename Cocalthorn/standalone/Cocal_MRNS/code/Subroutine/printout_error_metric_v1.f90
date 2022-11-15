subroutine printout_error_metric_v1(fnchar, iter_count,error_pot)
  implicit none
  real(8) :: error_pot
  integer :: iter_count
  character(len=4) :: fnchar
  write(6,'(a11,i4,a18,a4,a3,1p,e14.6)') 'Iteration= ', iter_count, &
  &      ', Error in metric ', fnchar, ' = ', error_pot
end subroutine printout_error_metric_v1
