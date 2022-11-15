subroutine interpf_ini(nrftmp,rftmp,flv)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrf
  use coordinate_grav_r_1D, only : rg
  implicit none
!
  real(8), intent(inout) :: rftmp(0:nnrg), flv(0:nnrg)
  integer, intent(in)    :: nrftmp
  real(8) :: fac, flvbk(0:nnrg)
  integer :: irf, irr, irf0
!
  flvbk(0:nrftmp) = flv(0:nrftmp)
  flv(1:nrftmp) = 0.0d0
!  flv(0) = flvbk(0)
!
  do irf = 1, nrf
    irf0 = 0
    do irr = 1, nrftmp
      if (rg(irf) >= rftmp(irr-1).and.rg(irf) < rftmp(irr)) then
        fac = (rg(irf) - rftmp(irr-1))/(rftmp(irr) - rftmp(irr-1))
        flv(irf) = fac*flvbk(irr) + (1.0d0 - fac)*flvbk(irr-1)
        irf0 = 1
        exit
      end if
    end do
    if (irf0 == 1) cycle
    if (rg(irf) >= rftmp(nrftmp)) flv(irf) = flvbk(nrftmp)
  end do
!
end subroutine interpf_ini
