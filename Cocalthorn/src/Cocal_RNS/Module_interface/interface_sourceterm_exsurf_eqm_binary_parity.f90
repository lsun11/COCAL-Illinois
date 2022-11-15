module interface_sourceterm_exsurf_eqm_binary_parity
  implicit none
  interface 
    subroutine sourceterm_exsurf_eqm_binary_parity(fnc,sou_ex,dsou_ex,par)
      real(8), pointer :: fnc(:,:,:), sou_ex(:,:), dsou_ex(:,:)
      real(8), intent(in) :: par
    end subroutine sourceterm_exsurf_eqm_binary_parity
  end interface
end module interface_sourceterm_exsurf_eqm_binary_parity
