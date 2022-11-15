subroutine calc_vector_phi_grav(sb)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use coordinate_grav_r, only : rg, hrg
  use trigonometry_grav_theta, only : sinthg, hsinthg
  use trigonometry_grav_phi, only : sinphig, cosphig, hsinphig, hcosphig
  use def_vector_phi, only : vec_phig, hvec_phig
  use def_binary_parameter, only : dis
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long) :: rrgg, sth, sphi, cphi, dis_bin
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
        sphi = sinphig(ipg)
        cphi = cosphig(ipg)
        vec_phig(irg,itg,ipg,1) =           - rrgg*sth*sphi
        vec_phig(irg,itg,ipg,2) = - dis_bin + rrgg*sth*cphi
        vec_phig(irg,itg,ipg,3) = 0.0d0
        if (irg.eq.0.or.itg.eq.0.or.ipg.eq.0) cycle
        rrgg = hrg(irg)
        sth  = hsinthg(itg)
        sphi = hsinphig(ipg)
        cphi = hcosphig(ipg)
        hvec_phig(irg,itg,ipg,1) =           - rrgg*sth*sphi
        hvec_phig(irg,itg,ipg,2) = - dis_bin + rrgg*sth*cphi
        hvec_phig(irg,itg,ipg,3) = 0.0d0
      end do
    end do
  end do
end subroutine calc_vector_phi_grav
