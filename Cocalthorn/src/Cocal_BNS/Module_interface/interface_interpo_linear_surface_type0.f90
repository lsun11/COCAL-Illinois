module interface_interpo_linear_surface_type0
  implicit none
  interface 
    subroutine interpo_linear_surface_type0(val,fnc,ir,it,ip)
      real(8) :: val
      real(8), pointer :: fnc(:,:,:)
      integer :: ir, it, ip
    end subroutine interpo_linear_surface_type0
  end interface
end module interface_interpo_linear_surface_type0
