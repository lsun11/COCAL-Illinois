module interface_error_matter
  implicit none
  interface 
    subroutine error_matter(pot,pot_bak,error,flag)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
    end subroutine error_matter
  end interface
end module interface_error_matter
