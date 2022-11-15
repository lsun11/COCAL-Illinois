subroutine IO_input_initial_vrot_NS_mpt(impt)
  use phys_constant, only  : long, pi
  use grid_parameter
  use coordinate_grav_r
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig
  use def_binary_parameter,    only : sepa, dis
  use def_metric
  use def_velocity_rot
  use def_matter_parameter
  implicit none
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, xycm, tcm, pcm, xcm,ycm,zcm
  real(long) :: rr, work_shift, pari
  integer    :: irf, itf, ipf, impt, npf_l, npf_r
  character(30) :: char1, char2, char3, char4, char5
!
!  write( 6,*)  sepa, dis

  confpow = 0.0d0
!

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
        rcm2 = (xa-0.5d0*sepa)**2 + ya**2 + za**2
        xycm2= (xa-0.5d0*sepa)**2 + ya**2
        rcm = dsqrt(rcm2)
        xycm= dsqrt(xycm2)
!       tcm = atan2(sqrt(xycm2),za)
!       pcm = dmod(2.0d0*pi+datan2(ya,xa-0.5d0*sepa),2.0d0*pi)
!        xcm = xa - 0.5d0*sepa
!        ycm = ya
!        zcm = za
        wxspf(irf,itf,ipf) = omespx*0.0d0 + omespy*(+za) + omespz*(-ya)
        wyspf(irf,itf,ipf) = omespx*(-za) + omespy*0.0d0 + omespz*(+xa)
        wzspf(irf,itf,ipf) = omespx*(+ya) + omespy*(-xa) + omespz*0.0d0
      end do
    end do
  end do
!
end subroutine IO_input_initial_vrot_NS_mpt
