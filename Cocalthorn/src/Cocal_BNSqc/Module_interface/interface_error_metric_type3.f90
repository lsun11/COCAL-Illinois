module interface_error_metric_type3
  implicit none
  interface 
    subroutine error_metric_type3(char_mp,pot,pot_bak,error,flag,ctype)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
      character(len=2), intent(in) :: ctype
      character(len=4), intent(in) :: char_mp
    end subroutine error_metric_type3
  end interface
end module interface_error_metric_type3
