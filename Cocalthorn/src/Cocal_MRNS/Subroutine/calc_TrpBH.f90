subroutine calc_TrpBH
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : rg, hrg
  use def_metric, only : psi
  use def_metric_pBH, only : R_trpBH,    hR_trpBH, psi_trpBH, hpsi_trpBH, &
  &                       alph_trpBH, halph_trpBH, bvr_trpBH, hbvr_trpBH
  use def_bh_parameter, only : mass_pBH
  implicit none
  integer    :: irg, itg, ipg
  real(long) :: rgc, rini, rval, faca, facb, MoverR
!
! --- Generate numerical values of analytic trunpet BH solution.
!
  faca = 27.0d0/16.0d0
  facb = 3.0d0*dsqrt(3.0d0)/4.0d0
  R_trpBH(0) = 3.0d0*mass_pBH/2.0d0
  psi_trpBH(0)  = 1.0d+20
  alph_trpBH(0) = 0.0d+00
  bvr_trpBH(0)  = 0.0d+00
  do irg = 1, nrg
    rini = R_trpBH(irg-1)
!
    rgc = rg(irg)
    call calc_eqsolver_TrpBH(rgc,rini,rval)
    MoverR = mass_pBH/rval
    R_trpBH(irg)   = rval
    psi_trpBH(irg) = dsqrt(rval/rgc)
    alph_trpBH(irg)= dsqrt(1.0d0 - 2.0d0*MoverR + faca*MoverR**4)
    bvr_trpBH(irg) = facb*MoverR**2*rgc/rval
!
    rgc = hrg(irg)
    call calc_eqsolver_TrpBH(rgc,rini,rval)
    MoverR = mass_pBH/rval
    hR_trpBH(irg)   = rval
    hpsi_trpBH(irg) = dsqrt(rval/rgc)
    halph_trpBH(irg)= dsqrt(1.0d0 - 2.0d0*MoverR + faca*MoverR**4)
    hbvr_trpBH(irg) = facb*MoverR**2*rgc/rval
  end do
!
!do ipg = 0, npg
!do itg = 0, ntg
!psi(0:nrg,itg,ipg) = psi_trpBH(0:nrg)
!end do
!end do
!do irg = 0, nrg
!write(6,*)rg(irg), R_trpBH(irg), psi_trpBH(irg)
!end do
!stop
end subroutine calc_TrpBH
