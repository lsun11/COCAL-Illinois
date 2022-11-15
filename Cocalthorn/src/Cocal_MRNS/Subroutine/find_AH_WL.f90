subroutine find_AH_WL
  use grid_parameter, only : iter_max
  implicit none
  integer :: iseq, iter_count, total_iteration
!
  call iteration_AHfinder_WL(iter_count)
  if (total_iteration.ge.iter_max) then
    write(6,*)' ** Solution did not converge **'
  end if
!
!  call calc_AHarea_AHfinder
  call calc_AHarea_WL

  call IO_output_AHfinder
  call IO_output_AHfinder_gnuplot
!
end subroutine find_AH_WL
