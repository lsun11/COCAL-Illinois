subroutine calc_potential_minus_rinv_excision
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use grid_parameter_binary_excision, only : 
  use coordinate_grav_r, only : rg
  use def_metric, only : psi
  use def_vector_x, only : vec_xg
  implicit none
  real(long) :: rrgg, small = 1.0d-14
  integer :: irg,itg,ipg
!
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        rrgg = sqrt(vec_xg(irg,itg,ipg,1)**2 + vec_xg(irg,itg,ipg,2)**2 &
        &         + vec_xg(irg,itg,ipg,3)**2 + small**2)
        if(rrgg >= 10.0d0) &
        & psi(irg,itg,ipg) = psi(irg,itg,ipg) + (2.0/(4.0*pi))/rrgg
      end do
    end do
  end do
!
!
end subroutine calc_potential_minus_rinv_excision
