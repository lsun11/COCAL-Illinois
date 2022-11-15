module interface_interpo_linear_type0_2Dsurf
  implicit none
  interface 
    subroutine interpo_linear_type0_2Dsurf(val,fnc,it,ip)
      real(8) :: val
      real(8), pointer :: fnc(:,:)
      integer :: it, ip
    end subroutine interpo_linear_type0_2Dsurf
  end interface
end module interface_interpo_linear_type0_2Dsurf
