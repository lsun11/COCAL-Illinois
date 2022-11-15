module interface_calc_dx_vec
  implicit none
  interface 
    subroutine calc_dx_vec(cline,ir,it,ip,dxv)
      character(len=2), intent(in) :: cline
      integer,          intent(in) :: ir,it,ip
      real(8), pointer             :: dxv(:)
    end subroutine calc_dx_vec
  end interface
end module interface_calc_dx_vec
