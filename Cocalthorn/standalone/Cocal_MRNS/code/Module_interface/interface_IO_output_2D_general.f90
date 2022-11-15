module interface_IO_output_2D_general
  implicit none
  interface 
    subroutine IO_output_2D_general(filename,fg,mg,pot,pl)
      real(8), pointer :: pot(:,:,:)
      character(len=2), intent(in) :: pl
      character(len=1), intent(in) :: fg,mg
      character(len=30), intent(in) :: filename
    end subroutine IO_output_2D_general
  end interface
end module interface_IO_output_2D_general
