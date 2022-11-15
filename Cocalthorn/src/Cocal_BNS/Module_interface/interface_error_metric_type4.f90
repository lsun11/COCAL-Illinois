module interface_error_metric_type4
  implicit none
  interface 
    subroutine error_metric_type4(pot,pot_bak,error,flag,ctype,cpot)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
      character(len=2), intent(in) :: ctype
      character(len=4), intent(in) :: cpot
    end subroutine error_metric_type4
  end interface
end module interface_error_metric_type4
