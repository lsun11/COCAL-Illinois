module interface_source_circ_surf_peos_spin
  implicit none
  interface 
    subroutine source_circ_surf_peos_spin(cline, ia, ib, souf)
      character(len=2), intent(in) :: cline
      integer,          intent(in) :: ia, ib
      real(8),          pointer    :: souf(:,:)
    end subroutine source_circ_surf_peos_spin
  end interface
end module interface_source_circ_surf_peos_spin
