subroutine copy_Aij_pBH_to_tfkij
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric,     only : psi, tfkij, tfkijkij
  use def_metric_excurve_grid, only : tfkij_grid, tfkijkij_grid
  use def_metric_pBH, only : aij_trpBH,      aijaij_trpBH, &
  &                          aij_trpBH_grid, aijaij_trpBH_grid
  use interface_interpo_linear_type0
  implicit none
  real(long) :: psim6, psim12, psigc
  integer :: irg, itg, ipg
!
! --- Compute components of extringic curvature 
! --- whose values are assigned on the mid points. 
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        psim6  = 1.0d0/psi(irg,itg,ipg)**6
        psim12 = 1.0d0/psi(irg,itg,ipg)**12
!
        tfkij_grid(irg,itg,ipg,1:3,1:3) = psim6 &
        &                               * aij_trpBH_grid(irg,itg,ipg,1:3,1:3) 
        tfkijkij_grid(irg,itg,ipg) = psim12*aijaij_trpBH_grid(irg,itg,ipg)
!
        if (irg.eq.0.or.itg.eq.0.or.ipg.eq.0) cycle
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        psim6  = 1.0d0/psigc**6
        psim12 = 1.0d0/psigc**12
        tfkij(irg,itg,ipg,1:3,1:3) = psim6*aij_trpBH(irg,itg,ipg,1:3,1:3)
        tfkijkij(irg,itg,ipg)  = psim12*aijaij_trpBH(irg,itg,ipg)
      end do
    end do
  end do
!
end subroutine copy_Aij_pBH_to_tfkij
