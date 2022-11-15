subroutine interpg_ini_iso(nrgtmp,rgtmp,grv)
!
  use phys_constant, only : nnrg
  use grid_parameter_1D, only : nrf
  use coordinate_grav_r_1D, only : rg
  implicit none
!
  real(8), intent(inout) :: rgtmp(0:nnrg), grv(0:nnrg)
  integer, intent(in)    :: nrgtmp
  real(8) :: x4(4), f4(4), fac, grvbk(0:nnrg)
  integer :: irg, irr, irg0
!
  if (nrgtmp>nnrg) then
    write(6,*) "Import_isotr_peos -> interpg_ini_iso"
    write(6,*) "Number of lines in psi, alph, emdg larger than maximum...exiting"
    stop
  end if

  grvbk(0:nrgtmp) = grv(0:nrgtmp)
  grv(1:nrgtmp) = 0.0d0
!  grv(0) = grvbk(0)
!
! For the isotropic solver irg=1,nrf. For the regular irg=1,nrg
  do irg = 1, nrf
    irg0 = 0
    do irr = 1, nrgtmp
      if (rg(irg) >= rgtmp(irr-1).and.rg(irg) < rgtmp(irr)) then
        fac = (rg(irg) - rgtmp(irr-1))/(rgtmp(irr) - rgtmp(irr-1))
        grv(irg) = fac*grvbk(irr) + (1.0d0 - fac)*grvbk(irr-1)
        irg0 = 1
        exit!go to 100
      end if
    end do
    if (irg0 == 1) cycle
    if (rg(irg) >= rgtmp(nrgtmp)) grv(irg) = grvbk(nrgtmp)
  end do
! 100  continue
!
end subroutine interpg_ini_iso
