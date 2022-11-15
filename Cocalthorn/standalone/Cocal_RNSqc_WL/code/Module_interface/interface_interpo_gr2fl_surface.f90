module interface_interpo_gr2fl_surface
  implicit none
  interface 
    subroutine interpo_gr2fl_surface(a,b)
      real(8), pointer :: a(:,:,:), b(:,:)
    end subroutine interpo_gr2fl_surface
  end interface
end module interface_interpo_gr2fl_surface
