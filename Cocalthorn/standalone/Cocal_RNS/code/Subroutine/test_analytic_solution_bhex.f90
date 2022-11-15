subroutine test_analytic_solution_bhex
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf, rgin
  use def_metric, only  :   bvxd, bvyd
  use coordinate_grav_r, only : rg
  use grid_points_binary_excision, only : rb
  implicit none
  integer     ::   irg,itg,ipg
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   mass 
!
  mass = 2.0d0*rgin
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        bvxd(irg,itg,ipg) = 1.0d0 + mass/(2.0d0*rg(irg)) &
      &                           + mass/(2.0d0*rb(irg,itg,ipg))
        bvyd(irg,itg,ipg) = 1.0d0 - mass/(2.0d0*rg(irg)) &
      &                           - mass/(2.0d0*rb(irg,itg,ipg))
      end do
    end do
  end do
!
end subroutine test_analytic_solution_bhex
