module interface_msec_store_f_vector
  implicit none
  interface 
    subroutine msec_store_f_vector(f_vector)
      real(8), pointer :: f_vector(:)
    end subroutine msec_store_f_vector
  end interface
end module interface_msec_store_f_vector
