subroutine calc_vector_x_matter(sb)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_matter, only  : rs
  use coordinate_grav_r, only : rg, hrg
  use trigonometry_grav_theta, only : sinthg, hsinthg, costhg, hcosthg
  use trigonometry_grav_phi, only : sinphig, cosphig, hsinphig, hcosphig
  use def_vector_x, only : vec_xf, hvec_xf
  use def_binary_parameter, only : dis
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long) :: rrgg, rrss, sth, cth, sphi, cphi, dis_bin
  integer :: irf, itf, ipf, sb
!
  dis_bin = dis
  if(sb.eq.1) dis_bin = 0.0d0   ! for a single rotating star
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        rrss = rs(itf,ipf)
        rrgg = rg(irf)
        sth  = sinthg(itf)
        cth  = costhg(itf)
        sphi = sinphig(ipf)
        cphi = cosphig(ipf)
        vec_xf(irf,itf,ipf,1) = - dis_bin + rrss*rrgg*sth*cphi
        vec_xf(irf,itf,ipf,2) =             rrss*rrgg*sth*sphi
        vec_xf(irf,itf,ipf,3) =             rrss*rrgg*cth
        if (irf.eq.0.or.itf.eq.0.or.ipf.eq.0) cycle
        call interpo_linear_type0_2Dsurf(rrss,rs,itf,ipf)
        rrgg = hrg(irf)
        sth  = hsinthg(itf)
        cth  = hcosthg(itf)
        sphi = hsinphig(ipf)
        cphi = hcosphig(ipf)
        hvec_xf(irf,itf,ipf,1) = - dis_bin + rrss*rrgg*sth*cphi
        hvec_xf(irf,itf,ipf,2) =             rrss*rrgg*sth*sphi
        hvec_xf(irf,itf,ipf,3) =             rrss*rrgg*cth
      end do
    end do
  end do
end subroutine calc_vector_x_matter
