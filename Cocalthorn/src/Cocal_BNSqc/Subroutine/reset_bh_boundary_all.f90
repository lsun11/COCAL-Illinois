subroutine reset_bh_boundary_all
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgin
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig
  use def_binary_parameter,    only : sepa
  use def_metric 
  implicit none
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, tcm, pcm
  integer    :: itg, ipg
!
  do ipg = 0, npg
    do itg = 0, ntg
      st = sinthg(itg)
      ct = costhg(itg)
      sp = sinphig(ipg)
      cp = cosphig(ipg)
! Coordinates wrt the center of black hole
      xa = rgin*st*cp
      ya = rgin*st*sp
      za = rgin*ct
! Coordinates wrt the CM
      rcm2 = (xa-0.5d0*sepa)**2 + ya**2 + za**2
      xycm2= (xa-0.5d0*sepa)**2 + ya**2

      rcm = sqrt(rcm2)
      tcm = atan2(sqrt(xycm2),za)
      pcm = dmod(2.0d0*pi+datan2(ya,xa-0.5d0*sepa),2.0d0*pi)
!
      psi(0,itg,ipg)  = 3.0d0
      alps(0,itg,ipg) = 1.0d0
      bvxd(0,itg,ipg) = -0.08d0*(-rcm*sin(tcm)*sin(pcm))
      bvyd(0,itg,ipg) = -0.08d0*(+rcm*sin(tcm)*cos(pcm))
      bvzd(0,itg,ipg) = 0.0d0
    end do
  end do
end subroutine reset_bh_boundary_all
