module interface_interpo_linear_type0_meridian
  implicit none
  interface 
    subroutine interpo_linear_type0_meridian(val,fnc,ir,it,ip)
      real(8) :: val
      real(8), pointer :: fnc(:,:,:)
      integer :: ir, it, ip
    end subroutine interpo_linear_type0_meridian
  end interface
end module interface_interpo_linear_type0_meridian
