module interface_hydrostatic_eq_qeos
  implicit none
  interface 
    subroutine hydrostatic_eq_qeos(rho)
      real(8), pointer :: rho(:,:,:)
    end subroutine hydrostatic_eq_qeos
  end interface
end module interface_hydrostatic_eq_qeos
