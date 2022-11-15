subroutine calc_vector_x_grav(sb)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use coordinate_grav_r, only : rg, hrg
  use trigonometry_grav_theta, only : sinthg, hsinthg, costhg, hcosthg
  use trigonometry_grav_phi, only : sinphig, cosphig, hsinphig, hcosphig
  use def_vector_x, only : vec_xg, hvec_xg
  use def_binary_parameter, only : dis
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long) :: rrgg, sth, cth, sphi, cphi, dis_bin
  integer :: irg, itg, ipg, sb
!
  dis_bin = dis
  if(sb.eq.1) dis_bin = 0.0d0   ! for a single rotating star
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        rrgg = rg(irg)
        sth  = sinthg(itg)
        cth  = costhg(itg)
        sphi = sinphig(ipg)
        cphi = cosphig(ipg)
        vec_xg(irg,itg,ipg,1) = - dis_bin + rrgg*sth*cphi
        vec_xg(irg,itg,ipg,2) =             rrgg*sth*sphi
        vec_xg(irg,itg,ipg,3) =             rrgg*cth
        if (irg.eq.0.or.itg.eq.0.or.ipg.eq.0) cycle
        rrgg = hrg(irg)
        sth  = hsinthg(itg)
        cth  = hcosthg(itg)
        sphi = hsinphig(ipg)
        cphi = hcosphig(ipg)
        hvec_xg(irg,itg,ipg,1) = - dis_bin + rrgg*sth*cphi
        hvec_xg(irg,itg,ipg,2) =             rrgg*sth*sphi
        hvec_xg(irg,itg,ipg,3) =             rrgg*cth
      end do
    end do
  end do
end subroutine calc_vector_x_grav
