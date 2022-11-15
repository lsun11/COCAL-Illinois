module interface_error_metric_type0
  implicit none
  interface 
    subroutine error_metric_type0(pot,pot_bak,error,flag,ctype)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
      character(len=2), intent(in) :: ctype
    end subroutine error_metric_type0
  end interface
end module interface_error_metric_type0
