module interface_error_metric
  implicit none
  interface 
    subroutine error_metric(pot,pot_bak,error,flag)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
    end subroutine error_metric
  end interface
end module interface_error_metric
