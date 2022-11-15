module interface_IO_output_1D_general
  implicit none
  interface 
    subroutine IO_output_1D_general(filename,fg,mg,pot,jr,jt,jp)
      real(8), pointer :: pot(:,:,:)
      integer, intent(in) :: jr, jt, jp
      character(len=1), intent(in) :: fg,mg
      character(len=30), intent(in) :: filename
    end subroutine IO_output_1D_general
  end interface
end module interface_IO_output_1D_general
