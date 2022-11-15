module interface_error_metric_type2_mpt
  implicit none
  interface 
    subroutine error_metric_type2_mpt(pot,pot_bak,error,flag,ctype,impt)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
      integer, intent(in)  :: impt
      character(len=2), intent(in) :: ctype
    end subroutine error_metric_type2_mpt
  end interface
end module interface_error_metric_type2_mpt
