module interface_interpo_fl2gr
  implicit none
  interface 
    subroutine interpo_fl2gr(a,b)
      real(8), pointer :: a(:,:,:), b(:,:,:)
    end subroutine interpo_fl2gr
  end interface
end module interface_interpo_fl2gr
