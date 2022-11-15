module interface_interpo_fl2gr_linear
  implicit none
  interface 
    subroutine interpo_fl2gr_linear(a,b)
      real(8), pointer :: a(:,:,:), b(:,:,:)
    end subroutine interpo_fl2gr_linear
  end interface
end module interface_interpo_fl2gr_linear
