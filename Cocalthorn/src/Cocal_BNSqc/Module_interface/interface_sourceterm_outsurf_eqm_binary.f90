module interface_sourceterm_outsurf_eqm_binary
  implicit none
  interface 
    subroutine sourceterm_outsurf_eqm_binary(fnc,sou_out,dsou_out)
      real(8), pointer :: fnc(:,:,:), sou_out(:,:), dsou_out(:,:)
    end subroutine sourceterm_outsurf_eqm_binary
  end interface
end module interface_sourceterm_outsurf_eqm_binary
