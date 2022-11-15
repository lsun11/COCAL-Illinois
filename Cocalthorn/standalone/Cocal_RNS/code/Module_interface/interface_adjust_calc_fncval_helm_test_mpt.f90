module interface_adjust_calc_fncval_helm_test_mpt
  implicit none
  interface 
    subroutine adjust_calc_fncval_helm_test_mpt(niq,msec_f)
      integer :: niq
      real(8) :: msec_f(niq)
    end subroutine adjust_calc_fncval_helm_test_mpt
  end interface
end module interface_adjust_calc_fncval_helm_test_mpt
