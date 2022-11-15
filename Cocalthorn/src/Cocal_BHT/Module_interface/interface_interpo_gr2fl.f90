module interface_interpo_gr2fl
  implicit none
  interface 
    subroutine interpo_gr2fl(a,b)
      real(8), pointer :: a(:,:,:), b(:,:,:)
    end subroutine interpo_gr2fl
  end interface
end module interface_interpo_gr2fl
