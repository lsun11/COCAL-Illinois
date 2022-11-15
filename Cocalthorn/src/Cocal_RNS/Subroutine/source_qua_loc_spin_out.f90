subroutine source_qua_loc_spin_out(soug,irs)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf
  use make_array_2d
  use def_metric, only :  psi, tfkij
  use coordinate_grav_r, only : hrg
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use def_vector_phi, only : hvec_phig
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: soug(:,:), psi_irs(:,:)
  real(long) :: ni, vphi_cm, Aij_surf, val, psi6, work(2,2)
  real(long) :: x(2),y(2), v, hvec_x(3)
  integer    :: irg, itg, ipg, ia, ib, irs, ii
  real(long) :: rrgg, sth, cth, sphi, cphi
!
  call alloc_array2d(psi_irs,0,ntg,0,npg)
!
  psi_irs(0:ntg,0:npg) = psi(irs,0:ntg,0:npg)
  rrgg = hrg(irs)

  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,psi_irs,itg,ipg)
      psi6 = val**6
      soug(itg,ipg)=0.0d0

      sth  = hsinthg(itg)
      cth  = hcosthg(itg)
      sphi = hsinphig(ipg)                                               
      cphi = hcosphig(ipg)

      hvec_x(1) = rrgg*sth*cphi
      hvec_x(2) = rrgg*sth*sphi
      hvec_x(3) = rrgg*cth

      do ib = 1, 3
        do ia = 1, 3
          Aij_surf = tfkij(irs,itg,ipg,ia,ib)
          ni = +hvec_x(ib)/rrgg
          vphi_cm = hvec_phig(irs,itg,ipg,ia)
          soug(itg,ipg) = soug(itg,ipg) + Aij_surf*vphi_cm*ni
        end do
      end do
      soug(itg,ipg) = soug(itg,ipg)*psi6
    end do
  end do
!
  deallocate(psi_irs)
end subroutine source_qua_loc_spin_out

