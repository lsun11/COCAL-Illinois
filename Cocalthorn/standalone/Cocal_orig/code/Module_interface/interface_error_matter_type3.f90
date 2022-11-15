module interface_error_matter_type3
  implicit none
  interface 
    subroutine error_matter_type3(pot,pot_bak,error,flag,cpot)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
      character(len=4), intent(in) :: cpot
    end subroutine error_matter_type3
  end interface
end module interface_error_matter_type3
