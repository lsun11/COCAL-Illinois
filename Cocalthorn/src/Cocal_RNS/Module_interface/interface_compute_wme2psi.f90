module interface_compute_wme2psi
  implicit none
  interface 
    subroutine compute_wme2psi(pot)
      real(8), pointer :: pot(:,:,:)
    end subroutine compute_wme2psi
  end interface
end module interface_compute_wme2psi
