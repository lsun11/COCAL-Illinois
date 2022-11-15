module interface_error_metric_type1
  implicit none
  interface 
    subroutine error_metric_type1(pot,pot_bak,error,ire,ite,ipe,flag,ctype)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag,ire,ite,ipe
      character(len=2), intent(in) :: ctype
    end subroutine error_metric_type1
  end interface
end module interface_error_metric_type1
