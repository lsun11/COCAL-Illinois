subroutine update_spin_BNS(iter_count)
  use phys_constant, only : long
  use grid_parameter
  use def_matter_parameter
  use def_velocity_rot
  use def_matter, only : rs
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : costhg, sinthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  implicit none
  real(long) :: st, ct, sp, cp, xa,ya,za, rr, rrss
  integer    :: irf, itf, ipf, iter_count

  if(iter_count==1) then
    write (6,*) "Update spin: omespx,omespy,omespz :", omespx, omespy, omespz 
    write (6,*) "--------------------------confpow :", confpow
  end if
  do irf = 0, nrf
    do ipf = 0, npf
      do itf = 0, ntf
        st = sinthg(itf)
        ct = costhg(itf)
        sp = sinphig(ipf)
        cp = cosphig(ipf)
        rr = rg(irf)
        rrss = rs(itf,ipf)
! 
        xa = rrss*rr*st*cp
        ya = rrss*rr*st*sp
        za = rrss*rr*ct
! 
        wxspf(irf,itf,ipf) = omespx*0.0d0 + omespy*(+za) + omespz*(-ya)
        wyspf(irf,itf,ipf) = omespx*(-za) + omespy*0.0d0 + omespz*(+xa)
        wzspf(irf,itf,ipf) = omespx*(+ya) + omespy*(-xa) + omespz*0.0d0
      end do
    end do
  end do
! 
end subroutine update_spin_BNS
