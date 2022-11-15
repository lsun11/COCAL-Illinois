subroutine calc_vector_irg(sb,irg)
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf, nrg, ntg, npg
  use coordinate_grav_r, only : rg, hrg
  use trigonometry_grav_theta, only : sinthg, hsinthg, costhg, hcosthg
  use trigonometry_grav_phi, only : sinphig, cosphig, hsinphig, hcosphig
  use def_binary_parameter, only : dis
!  use def_bh_parameter, only : th_spin_bh, phi_spin_bh
  use def_vector_irg
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long) :: rrgg, sth, cth, sphi, cphi, dis_bin
!  real(long) :: spin_matrix(3,3),  spin_axis(3), spin_z(3)
  integer :: irg, itg, ipg, sb, ii
!
  dis_bin = dis
  if(sb.eq.1) dis_bin = 0.0d0   ! for a single bh
!
!  spin_matrix(1,1) =  dcos(th_spin_bh)*dcos(phi_spin_bh)
!  spin_matrix(1,2) =                  -dsin(phi_spin_bh)
!  spin_matrix(1,3) =  dsin(th_spin_bh)*dcos(phi_spin_bh)
!  spin_matrix(2,1) =  dcos(th_spin_bh)*dsin(phi_spin_bh)
!  spin_matrix(2,2) =                   dcos(phi_spin_bh)
!  spin_matrix(2,3) =  dsin(th_spin_bh)*dsin(phi_spin_bh)
!  spin_matrix(3,1) = -dsin(th_spin_bh)
!  spin_matrix(3,2) =  0.0d0
!  spin_matrix(3,3) =  dcos(th_spin_bh)
!  spin_z(1) = 0.0d0
!  spin_z(2) = 0.0d0
!  spin_z(3) = 1.0d0
!  do ii = 1, 3
!    spin_axis(ii) = spin_matrix(ii,1)*spin_z(1) &
!    &             + spin_matrix(ii,2)*spin_z(2) &
!    &             + spin_matrix(ii,3)*spin_z(3)
!  end do
!
  do ipg = 0, npg
    do itg = 0, ntg
!      irg = 0
      rrgg = rg(irg)
      sth  = sinthg(itg)
      cth  = costhg(itg)
      sphi = sinphig(ipg)
      cphi = cosphig(ipg)
      vec_irg_cm_xg(itg,ipg,1) = - dis_bin + rrgg*sth*cphi
      vec_irg_cm_xg(itg,ipg,2) =             rrgg*sth*sphi
      vec_irg_cm_xg(itg,ipg,3) =             rrgg*cth
      vec_irg_cm_phig(itg,ipg,1) =           - rrgg*sth*sphi
      vec_irg_cm_phig(itg,ipg,2) = - dis_bin + rrgg*sth*cphi
      vec_irg_cm_phig(itg,ipg,3) = 0.0d0
      vec_irg_cbh_xg(itg,ipg,1) = rrgg*sth*cphi
      vec_irg_cbh_xg(itg,ipg,2) = rrgg*sth*sphi
      vec_irg_cbh_xg(itg,ipg,3) = rrgg*cth
      vec_irg_cbh_phig(itg,ipg,1) = - rrgg*sth*sphi
      vec_irg_cbh_phig(itg,ipg,2) =   rrgg*sth*cphi
      vec_irg_cbh_phig(itg,ipg,3) = 0.0d0
!      vec_irg_cbh_spin(itg,ipg,1) = spin_axis(2)*vec_irg_cbh_xg(itg,ipg,3)&
!      &                          - spin_axis(3)*vec_irg_cbh_xg(itg,ipg,2)
!      vec_irg_cbh_spin(itg,ipg,2) = spin_axis(3)*vec_irg_cbh_xg(itg,ipg,1)&
!      &                          - spin_axis(1)*vec_irg_cbh_xg(itg,ipg,3)
!      vec_irg_cbh_spin(itg,ipg,3) = spin_axis(1)*vec_irg_cbh_xg(itg,ipg,2)&
!      &                          - spin_axis(2)*vec_irg_cbh_xg(itg,ipg,1)
      if (itg.eq.0.or.ipg.eq.0) cycle
      rrgg = rg(irg) ! irg surface
      sth  = hsinthg(itg)
      cth  = hcosthg(itg)
      sphi = hsinphig(ipg)
      cphi = hcosphig(ipg)
      hvec_irg_cm_xg(itg,ipg,1) = - dis_bin + rrgg*sth*cphi
      hvec_irg_cm_xg(itg,ipg,2) =             rrgg*sth*sphi
      hvec_irg_cm_xg(itg,ipg,3) =             rrgg*cth
      hvec_irg_cm_phig(itg,ipg,1) =           - rrgg*sth*sphi
      hvec_irg_cm_phig(itg,ipg,2) = - dis_bin + rrgg*sth*cphi
      hvec_irg_cm_phig(itg,ipg,3) = 0.0d0
      hvec_irg_cbh_xg(itg,ipg,1) = rrgg*sth*cphi
      hvec_irg_cbh_xg(itg,ipg,2) = rrgg*sth*sphi
      hvec_irg_cbh_xg(itg,ipg,3) = rrgg*cth
      hvec_irg_cbh_phig(itg,ipg,1) = - rrgg*sth*sphi
      hvec_irg_cbh_phig(itg,ipg,2) =   rrgg*sth*cphi
      hvec_irg_cbh_phig(itg,ipg,3) = 0.0d0
!      hvec_irg_cbh_spin(itg,ipg,1) = spin_axis(2)*hvec_irg_cbh_xg(itg,ipg,3)&
!      &                           - spin_axis(3)*hvec_irg_cbh_xg(itg,ipg,2)
!      hvec_irg_cbh_spin(itg,ipg,2) = spin_axis(3)*hvec_irg_cbh_xg(itg,ipg,1)&
!      &                           - spin_axis(1)*hvec_irg_cbh_xg(itg,ipg,3)
!      hvec_irg_cbh_spin(itg,ipg,3) = spin_axis(1)*hvec_irg_cbh_xg(itg,ipg,2)&
!      &                           - spin_axis(2)*hvec_irg_cbh_xg(itg,ipg,1)
    end do
  end do
end subroutine calc_vector_irg
