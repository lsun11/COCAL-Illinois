module interface_update_matter
  implicit none
  interface 
    subroutine update_matter(potf,mtfield,convf)
      real(8), pointer :: potf(:,:,:)
      real(8), pointer :: mtfield(:,:,:)
      real(8), intent(in) :: convf
    end subroutine update_matter
  end interface
end module interface_update_matter
