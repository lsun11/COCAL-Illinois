subroutine source_ang_mom_exc(sous)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use make_array_2d
  use def_metric, only :  psi, tfkij
  use def_metric_excurve_grid, only : tfkij_grid
  use coordinate_grav_r, only : rg
  use def_vector_irg, only : hvec_irg_cm_phig, hvec_irg_cbh_xg
  use interface_interpo_linear_type0_2Dsurf
  use grid_parameter_binary_excision, only: ex_nrg
  implicit none
  real(long), external :: lagint_2nd
  real(long), pointer :: sous(:,:), psi_exc(:,:)
  real(long) :: ni, vphi_bh, vphi_cm, Aij_surf, val, psi6, work(2,2)
  real(long) :: x(2),y(2), v
  integer    :: irg, itg, ipg, ia, ib
!
  call alloc_array2d(psi_exc,0,ntg,0,npg)
  call calc_vector_irg(2,ex_nrg)
!
!  do ia = 1,3
!    do ib = 1,3
!      do ipg = 0,npg
!        do itg = 0,ntg
!          x(1) = rg(1)
!          x(2) = rg(2)
!          y(1) = tfkij_grid(1,itg,ipg,ia,ib)
!          y(2) = tfkij_grid(2,itg,ipg,ia,ib)
!          v = rg(0)      
!          tfkij_grid(0,itg,ipg,ia,ib) = lagint_2nd(x,y,v)
!        end do
!      end do
!    end do
!  end do

  psi_exc(0:ntg,0:npg) = psi(ex_nrg,0:ntg,0:npg)
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,psi_exc,itg,ipg)
      psi6 = val**6
      sous(itg,ipg)=0.0d0
      do ib = 1, 3
        do ia = 1, 3
          work(1:2,1:2) = tfkij_grid(0,itg-1:itg,ipg-1:ipg,ia,ib)
          call interpo_linear1p_type0_2Dsurf(Aij_surf,work)          
!         Aij_surf = tfkij(1,itg,ipg,ia,ib)
!
          ni = -hvec_irg_cbh_xg(itg,ipg,ib)/rg(ex_nrg)
          vphi_cm = hvec_irg_cm_phig(itg,ipg,ia)
!
          sous(itg,ipg) = sous(itg,ipg) + Aij_surf*vphi_cm*ni
        end do
      end do
      sous(itg,ipg) = sous(itg,ipg)*psi6
    end do
  end do
!
  deallocate(psi_exc)
end subroutine source_ang_mom_exc

