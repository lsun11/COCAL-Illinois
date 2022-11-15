module interface_msec_copy_to_iter_quantities
  implicit none
  interface 
    subroutine msec_copy_to_iter_quantities(x_vector)
      real(8), pointer :: x_vector(:)
    end subroutine msec_copy_to_iter_quantities
  end interface
end module interface_msec_copy_to_iter_quantities
