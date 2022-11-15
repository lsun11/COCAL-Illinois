module interface_adjust_copy_ome_cm_ratio_to_mpt
  implicit none
  interface 
    subroutine adjust_copy_ome_cm_ratio_to_mpt(niq,msec_x)
      integer :: niq
      real(8) :: msec_x(niq)
    end subroutine adjust_copy_ome_cm_ratio_to_mpt
  end interface
end module interface_adjust_copy_ome_cm_ratio_to_mpt
