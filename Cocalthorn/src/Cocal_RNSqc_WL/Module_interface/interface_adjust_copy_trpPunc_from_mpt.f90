module interface_adjust_copy_trpPunc_from_mpt
  implicit none
  interface 
    subroutine adjust_copy_trpPunc_from_mpt(niq,msec_x)
      integer :: niq
      real(8) :: msec_x(niq)
    end subroutine adjust_copy_trpPunc_from_mpt
  end interface
end module interface_adjust_copy_trpPunc_from_mpt
