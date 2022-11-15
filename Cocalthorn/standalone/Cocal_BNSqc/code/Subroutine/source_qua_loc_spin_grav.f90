subroutine source_qua_loc_spin_grav(soug,irs)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf
  use make_array_2d
  use def_metric, only :  psi, tfkij
  use coordinate_grav_r, only : hrg
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use def_vector_phi, only : hvec_phig
  use interface_interpo_linear_type0 
  use make_array_2d
  implicit none
  real(long), pointer :: soug(:,:)
  real(long) :: ni, vphi_cm, Aij_surf, psi6
  real(long) :: x(2),y(2), val, hvec_x(3), phi_ns(3)
  integer    :: irg, itg, ipg, ia, ib, irs, ic
  real(long) :: sth, cth, sphi, cphi
!
!  irs is a midpoint         
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0(val,psi,irs,itg,ipg)
      psi6 = val**6
      soug(itg,ipg)=0.0d0

      sth  = hsinthg(itg)
      cth  = hcosthg(itg)
      sphi = hsinphig(ipg)                                               
      cphi = hcosphig(ipg)

      hvec_x(1) = sth*cphi
      hvec_x(2) = sth*sphi
      hvec_x(3) = cth

      phi_ns(1) = - hrg(irs)*sth*sphi
      phi_ns(2) =   hrg(irs)*sth*cphi
      phi_ns(3) =   0.0d0

      do ib = 1, 3
        do ia = 1, 3
          Aij_surf = tfkij(irs,itg,ipg,ia,ib)
!          vphi_cm = hvec_phig(irs,itg,ipg,ib)
          vphi_cm = phi_ns(ib)
          ni = +hvec_x(ia)
          soug(itg,ipg) = soug(itg,ipg) + Aij_surf*vphi_cm*ni
        end do
      end do
      soug(itg,ipg) = soug(itg,ipg)*psi6
    end do
  end do
!
end subroutine source_qua_loc_spin_grav

