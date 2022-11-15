subroutine test_analytic_BHNS_solution_mpt(impt)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   nrg, ntg, npg, nrf
  use def_metric, only  :   bvxd
  use coordinate_grav_r, only : rg
  use grid_points_binary_excision, only : rb
  implicit none
  integer     ::   irg, itg, ipg, impt
  real(long)  ::   zfac, small = 1.0d-15
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
      if(impt.eq.1) then
        bvxd(irg,itg,ipg) = -(8.0d0*pi/3.0d0)*1.0d0/rg(irg)  &
        &                     -(4.0d0*pi/3.0d0)*1.0d0/rb(irg,itg,ipg)
      end if
      if(impt.eq.2) then
        if(irg.lt.nrf) &
        & bvxd(irg,itg,ipg) = - 2.0d0*pi/3.0d0*(3.0d0-rg(irg)**2) &
        &                     -(8.0d0*pi/3.0d0)*1.0d0/rb(irg,itg,ipg)
        if(irg.ge.nrf) &
        & bvxd(irg,itg,ipg) = -(4.0d0*pi/3.0d0)*1.0d0/rg(irg)  &
        &                     -(8.0d0*pi/3.0d0)*1.0d0/rb(irg,itg,ipg)
      end if
!!      if(irg.eq.0) &
!!      & bvxd(irg,itg,ipg) = - 4.0d0/pi &
!!      &                     - 4.0d0/pi - 4.0d0/(pi*rb(irg,itg,ipg))
!!      if(irg.lt.nrf.and.irg.ne.0) &
!!      & bvxd(irg,itg,ipg) = - 4.0d0*sin(pi*rg(irg)) &
!!      &                      /(pi**2*rg(irg)) &
!!      &                     - 4.0d0/pi - 4.0d0/(pi*rb(irg,itg,ipg))
!!      if(irg.ge.nrf) &
!!      & bvxd(irg,itg,ipg) = - 4.0d0/(pi*rg(irg))  &
!!      &                     - 4.0d0/(pi*rb(irg,itg,ipg))
      end do
    end do
  end do
!
end subroutine test_analytic_BHNS_solution_mpt
