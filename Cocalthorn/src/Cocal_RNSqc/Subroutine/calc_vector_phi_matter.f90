subroutine calc_vector_phi_matter(sb)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use def_matter, only  : rs
  use coordinate_grav_r, only : rg, hrg
  use trigonometry_grav_theta, only : sinthg, hsinthg
  use trigonometry_grav_phi, only : sinphig, cosphig, hsinphig, hcosphig
  use def_vector_phi, only : vec_phif, hvec_phif, hvec_phif_surface
  use def_binary_parameter, only : dis
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long) :: rrgg, rrss, sth, sphi, cphi, dis_bin
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
        sphi = sinphig(ipf)
        cphi = cosphig(ipf)
        vec_phif(irf,itf,ipf,1) =           - rrss*rrgg*sth*sphi
        vec_phif(irf,itf,ipf,2) = - dis_bin + rrss*rrgg*sth*cphi
        vec_phif(irf,itf,ipf,3) = 0.0d0
        if (irf.eq.0.or.itf.eq.0.or.ipf.eq.0) cycle
        call interpo_linear_type0_2Dsurf(rrss,rs,itf,ipf)
        rrgg = hrg(irf)
        sth  = hsinthg(itf)
        sphi = hsinphig(ipf)
        cphi = hcosphig(ipf)
        hvec_phif(irf,itf,ipf,1) =           - rrss*rrgg*sth*sphi
        hvec_phif(irf,itf,ipf,2) = - dis_bin + rrss*rrgg*sth*cphi
        hvec_phif(irf,itf,ipf,3) = 0.0d0
        rrgg = rg(nrf)  ! at surface
        hvec_phif_surface(itf,ipf,1) =           - rrss*rrgg*sth*sphi
        hvec_phif_surface(itf,ipf,2) = - dis_bin + rrss*rrgg*sth*cphi
        hvec_phif_surface(itf,ipf,3) = 0.0d0
      end do
    end do
  end do
end subroutine calc_vector_phi_matter
