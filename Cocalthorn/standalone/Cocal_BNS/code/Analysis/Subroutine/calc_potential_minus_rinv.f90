subroutine calc_potential_minus_rinv
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg
  use def_metric, only : psi
  implicit none
  integer :: irg,itg,ipg
!
!
  do ipg = 0, npg
  do itg = 0, ntg
  do irg = 0, nrg
    if(rg(irg) >= 10.0d0) &
   & psi(irg,itg,ipg) = psi(irg,itg,ipg) + (2.0/(4.0*pi))/rg(irg)
  end do
  end do
  end do
!
!
end subroutine calc_potential_minus_rinv
