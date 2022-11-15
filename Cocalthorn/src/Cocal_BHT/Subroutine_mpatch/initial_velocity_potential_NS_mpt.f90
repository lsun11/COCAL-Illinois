subroutine initial_velocity_potential_NS_mpt(impt)
  use phys_constant, only  : long, pi
  use grid_parameter
  use coordinate_grav_r
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig
  use def_binary_parameter,    only : sepa, dis
  use def_matter_parameter, only  : ome
  use def_metric
  use def_matter
  use def_velocity_rot
  use def_velocity_potential
  implicit none
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, xycm, tcm, pcm, xcm,ycm,zcm
  real(long) :: rr, work_shift, pari, qq, hh, pre, rho0, ene
  integer    :: irf, itf, ipf, impt, npf_l, npf_r
  character(30) :: char1, char2, char3, char4, char5
!
!  write( 6,*)  sepa, dis

  if (impt.eq.2) then
    work_shift = dis;    pari = -1.0;
  else
    work_shift = -dis;   pari = 1.0;
  end if

  do irf = 0, nrf
    do ipf = 0, npf
      do itf = 0, ntf
        st = sinthg(itf)
        ct = costhg(itf)
        sp = sinphig(ipf)
        cp = cosphig(ipf)
        rr = rg(irf)
! Coordinates wrt the center of black hole
        xa = rr*st*cp
        ya = rr*st*sp
        za = rr*ct
! Coordinates wrt the CM
!        rcm2 = (xa-0.5d0*sepa)**2 + ya**2 + za**2
!        xycm2= (xa-0.5d0*sepa)**2 + ya**2
!        rcm = dsqrt(rcm2)
!        xycm= dsqrt(xycm2)
!       tcm = atan2(sqrt(xycm2),za)
!       pcm = dmod(2.0d0*pi+datan2(ya,xa-0.5d0*sepa),2.0d0*pi)
        xcm = (-dis + xa)*pari
        ycm = ya*pari
        zcm = 0.0d0
        qq = emd(irf,itf,ipf)
        call peos_q2hprho(qq, hh, pre, rho0, ene)

        vepxf(irf,itf,ipf) = 0.0d0
        vepyf(irf,itf,ipf) = ome*(-dis)
        vepzf(irf,itf,ipf) = 0.0d0

        vep(irf,itf,ipf) = vepyf(irf,itf,ipf)*ya
      end do
    end do
  end do
!
end subroutine initial_velocity_potential_NS_mpt
