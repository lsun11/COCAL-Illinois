module interface_sourceterm_exsurf_eqm_binary
  implicit none
  interface 
    subroutine sourceterm_exsurf_eqm_binary(fnc,sou_ex,dsou_ex)
      real(8), pointer :: fnc(:,:,:), sou_ex(:,:), dsou_ex(:,:)
    end subroutine sourceterm_exsurf_eqm_binary
  end interface
end module interface_sourceterm_exsurf_eqm_binary
