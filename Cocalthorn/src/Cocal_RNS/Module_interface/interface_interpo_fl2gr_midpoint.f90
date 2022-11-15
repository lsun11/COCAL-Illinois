module interface_interpo_fl2gr_midpoint
  implicit none
  interface 
    subroutine interpo_fl2gr_midpoint(a,b)
      real(8), pointer :: a(:,:,:), b(:,:,:)
    end subroutine interpo_fl2gr_midpoint
  end interface
end module interface_interpo_fl2gr_midpoint
