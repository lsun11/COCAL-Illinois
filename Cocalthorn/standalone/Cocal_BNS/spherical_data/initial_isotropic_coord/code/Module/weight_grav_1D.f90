module weight_grav_1D
  use phys_constant, only : nnrg
  implicit none
  real(8) :: wgdrg(nnrg), w4drg(0:nnrg)	!def_weigg
  real(8) :: wgdr(0:nnrg), w4dr(0:nnrg)	!def_weigh
  real(8) :: drf(0:nnrg)
contains
subroutine weight_calc_grav_1D
!
  use grid_parameter_1D, only : nrf, nrg
  use coordinate_grav_r_1D, only : hrg, rg, drg
  implicit none
  real(8) :: san1, san2, san4
  integer :: ir
!
! --- weight for GR integration.
  san1 = 1.0d0/3.0d0
  san4 = 4.0d0/3.0d0
  san2 = 2.0d0/3.0d0
!
  wgdrg(1:nrg) = hrg(1:nrg)**2*drg(1:nrg)
  w4drg(1:nrg) = wgdrg(1:nrg)
!
! --- weight for fluid integration.
!
  ir = 0
  wgdr(ir) = 0.5d0*rg(ir)**2*drg(ir+1)
  ir = nrf
  wgdr(ir) = 0.5d0*rg(ir)**2*drg(ir)
  wgdr(1:nrf-1) = 0.5d0*rg(1:nrf-1)**2*drg(1:nrf-1) &
  &             + 0.5d0*rg(1:nrf-1)**2*drg(2:nrf)
!
  do ir = 1, nrf, 4
    w4dr(ir  ) = 13.0d0/12.0d0*wgdrg(ir)
    w4dr(ir+1) = 11.0d0/12.0d0*wgdrg(ir+1)
    w4dr(ir+2) = 11.0d0/12.0d0*wgdrg(ir+2)
    w4dr(ir+3) = 13.0d0/12.0d0*wgdrg(ir+3)
  end do
!
  drf(1:nrf) = drg(1:nrf)
!
end subroutine weight_calc_grav_1D
end module weight_grav_1D
