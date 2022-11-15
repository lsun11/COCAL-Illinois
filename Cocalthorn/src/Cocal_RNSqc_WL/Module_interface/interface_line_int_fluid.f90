module interface_line_int_fluid
  implicit none
  interface 
    subroutine line_int_fluid(cline, ia, ib, souf, line_int)
      character(len=2), intent(in) :: cline
      integer,          intent(in) :: ia, ib
      real(8), pointer             :: souf(:,:)
      real(8), intent(out)         :: line_int
    end subroutine line_int_fluid
  end interface
end module interface_line_int_fluid
