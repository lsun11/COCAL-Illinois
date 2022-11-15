module interface_source_circ_line_peos_irrot
  implicit none
  interface 
    subroutine source_circ_line_peos_irrot(cline, ia, ib, souf)
      character(len=2), intent(in) :: cline
      integer,          intent(in) :: ia, ib
      real(8),          pointer    :: souf(:,:)
    end subroutine source_circ_line_peos_irrot
  end interface
end module interface_source_circ_line_peos_irrot
