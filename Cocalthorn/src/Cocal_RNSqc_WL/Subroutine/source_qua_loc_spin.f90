subroutine source_qua_loc_spin(sous)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use make_array_2d
  use def_metric, only :  psi, tfkij
  use def_metric_excurve_grid, only : tfkij_grid
  use coordinate_grav_r, only : rg
  use def_vector_bh, only : hvec_bh_cbh_phig, hvec_bh_cbh_xg
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sous(:,:), psi_bh(:,:)
  real(long) :: ni, vphi_bh, vphi_cm, Aij_surf, val, psi6, work(2,2)
  real(long) :: x(2),y(2), v
  integer    :: irg, itg, ipg, ia, ib
!
  call alloc_array2d(psi_bh,0,ntg,0,npg)
  call calc_vector_bh(2)
!
  psi_bh(0:ntg,0:npg) = psi(0,0:ntg,0:npg)
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,psi_bh,itg,ipg)
      psi6 = val**6
      sous(itg,ipg)=0.0d0
      do ib = 1, 3
        do ia = 1, 3
          work(1:2,1:2) = tfkij_grid(0,itg-1:itg,ipg-1:ipg,ia,ib)
          call interpo_linear1p_type0_2Dsurf(Aij_surf,work)          
!
          ni = -hvec_bh_cbh_xg(itg,ipg,ib)/rg(0)
          vphi_cm = hvec_bh_cbh_phig(itg,ipg,ia)
!
          sous(itg,ipg) = sous(itg,ipg) + Aij_surf*vphi_cm*ni
        end do
      end do
      sous(itg,ipg) = sous(itg,ipg)*psi6
    end do
  end do
!
  deallocate(psi_bh)
end subroutine source_qua_loc_spin

