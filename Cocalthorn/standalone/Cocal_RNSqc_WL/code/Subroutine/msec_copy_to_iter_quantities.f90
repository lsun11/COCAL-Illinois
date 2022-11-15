subroutine msec_copy_to_iter_quantities(x_vector)
  use phys_constant, only : long
  use def_bh_parameter
  use def_quantities
  use def_iter_quantities
  implicit none
  real(long), pointer  :: x_vector(:)
  integer :: i, niq, istep
!
  i=0
  i=i+1; ome_bh  = x_vector(i) 
  i=i+1; spin_bh = x_vector(i) 
!
end subroutine msec_copy_to_iter_quantities
