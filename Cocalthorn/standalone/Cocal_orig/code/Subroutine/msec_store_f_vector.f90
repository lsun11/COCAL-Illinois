subroutine msec_store_f_vector(f_vector)
  use phys_constant, only : long
  use def_bh_parameter
  use def_quantities
  use def_iter_quantities
  implicit none
  real(long), pointer :: f_vector(:)
  integer :: i
!
  i=0
  i=i+1; f_vector(i) = (admmass - komarmass)/admmass
  i=i+1; f_vector(i) = qua_loc_spin/irredmass/irredmass - 0.0d0
!
end subroutine msec_store_f_vector
