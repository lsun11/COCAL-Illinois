subroutine gravmid_sph(sou,pot)
!
  use phys_constant, only : nnrg
  use weight_grav_1D, only : wgdrg
  use grid_parameter_1D, only : nrg
  use radial_perm_fn_grav_1D, only : hfsn
  implicit none
!
  real(8), intent(inout) :: sou(0:nnrg), pot(0:nnrg)
  integer :: ir, irr
!
! --- Gravitational potential is caluculated. ---
!
!     L[p] = S --> p = 1/4pi Int[ S/r dV]
!              
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!    hfsn(irr,ir) -->  for r'< r
!                      for r < r'
!
  pot(0:nrg) = 0.0d0
  do ir = 0, nrg
    do irr = 1, nrg
!      wei = wgdrg(irr)
      pot(ir) = pot(ir) + sou(irr)*hfsn(irr,ir)*wgdrg(irr)
    end do
  end do
!
  pot(0:nrg) = - pot(0:nrg)
!
end subroutine gravmid_sph
