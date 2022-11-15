subroutine compute_alps2alph(pot,psi)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  implicit none
  real(long), pointer    :: pot(:,:,:), psi(:,:,:)
  integer :: irg, itg, ipg
  do irg = 0, nrg
    do itg = 0, ntg
      do ipg = 0, npg
        pot(irg, itg, ipg) = pot(irg, itg, ipg)/psi(irg, itg, ipg)
      end do
    end do
  end do
end subroutine compute_alps2alph
