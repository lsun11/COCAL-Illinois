subroutine test_analytic_solution_1bh_psialph
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgin
  use def_metric, only     : bvxd, bvyd
  use coordinate_grav_r, only : rg
  implicit none
  integer    :: irg, itg, ipg
  real(long) :: zfac, small = 1.0d-15
  real(long) :: mass 
!
  mass = 2.0d0*rgin
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        bvxd(irg,itg,ipg) = 1.0d0 + mass/(2.0d0*rg(irg))
        bvyd(irg,itg,ipg) =(1.0d0 - mass/(2.0d0*rg(irg)))/bvxd(irg,itg,ipg)
      end do
    end do
  end do
!
end subroutine test_analytic_solution_1bh_psialph
