module interface_compute_alps2alph
  implicit none
  interface 
    subroutine compute_alps2alph(pot,psi)
      real(8), pointer     :: pot(:,:,:), psi(:,:,:)
    end subroutine compute_alps2alph
  end interface
end module interface_compute_alps2alph
