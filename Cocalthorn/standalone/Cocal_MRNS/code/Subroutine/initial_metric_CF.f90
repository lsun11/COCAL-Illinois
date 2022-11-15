subroutine initial_metric_CF
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgin, ntgxy, npgxzm
  use coordinate_grav_r
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi,   only : sinphig, cosphig
  use def_binary_parameter,    only : sepa
  use def_metric
  implicit none
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, tcm, pcm, xcm,ycm,zcm
  real(long) :: rr
  integer    :: irg, itg, ipg
!
  do irg = 0, nrg
    do ipg = 0, npg
      do itg = 0, ntg
        st = sinthg(itg)
        ct = costhg(itg)
        sp = sinphig(ipg)
        cp = cosphig(ipg)
        rr = rg(irg)
! Coordinates wrt the center of black hole                                                                              
        xa = rr*st*cp
        ya = rr*st*sp
        za = rr*ct
! Coordinates wrt the CM                                                                                                 
!       rcm2 = (xa-0.5d0*sepa)**2 + ya**2 + za**2                                                                        
!       xycm2= (xa-0.5d0*sepa)**2 + ya**2                                                                                 
!       rcm = sqrt(rcm2)                                                                                                  
!       tcm = atan2(sqrt(xycm2),za)                                                                                       
!       pcm = dmod(2.0d0*pi+datan2(ya,xa-0.5d0*sepa),2.0d0*pi)                                                            
        xcm = xa - 0.5d0*sepa
        ycm = ya
        zcm = za
        psi(irg,itg,ipg)  = 1.0d0
        alph(irg,itg,ipg) = 1.0d0
        alps(irg,itg,ipg) = 1.0d0
        bvxd(irg,itg,ipg) = 0.0d0
        bvyd(irg,itg,ipg) = 0.0d0
        bvzd(irg,itg,ipg) = 0.0d0
!        bvxd(irg,itg,ipg) = st*sp
!        bvyd(irg,itg,ipg) = -st*cp
!        bvzd(irg,itg,ipg) = ct
      end do
    end do
  end do
!
end subroutine initial_metric_CF
