module interface_error_metric_type3_mpt
  implicit none
  interface 
    subroutine error_metric_type3_mpt(pot,pot_bak,error,flag,ctype,impt,cpot)
      real(8), pointer     :: pot(:,:,:), pot_bak(:,:,:)
      real(8), intent(out) :: error
      integer, intent(out) :: flag
      integer, intent(in)  :: impt
      character(len=2), intent(in) :: ctype
      character(len=4), intent(in) :: cpot
    end subroutine error_metric_type3_mpt
  end interface
end module interface_error_metric_type3_mpt
