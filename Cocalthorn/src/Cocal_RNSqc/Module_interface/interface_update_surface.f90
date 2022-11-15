module interface_update_surface
  implicit none
  interface 
    subroutine update_surface(potrs,rsnew,convf)
      real(8), pointer :: potrs(:,:)
      real(8), pointer :: rsnew(:,:)
      real(8), intent(in) :: convf
    end subroutine update_surface
  end interface
end module interface_update_surface
