subroutine source_ang_mom_smarr(sous)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use make_array_2d
  use def_metric
  use def_metric_excurve_grid, only : tfkij_grid
  use def_bh_parameter, only : ome_bh
  use coordinate_grav_r, only : rg
  use interface_interpo_linear_type0_2Dsurf
  use def_vector_bh, only : hvec_bh_cbh_xg, hvec_bh_cm_phig
  use interface_sourceterm_surface_int
  implicit none
  real(long), pointer :: sous(:,:), sous1(:,:), sous2(:,:), psi_bh(:,:)
  real(long), pointer :: bvxd_bh(:,:), bvyd_bh(:,:), bvzd_bh(:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  real(long) :: ni, Aij_surf, val, psi2, psi6, work(2,2), bv(3), beta, vphi_cm
  integer    :: irg, itg, ipg, ia, ib
!
  call alloc_array2d(psi_bh,0,ntg,0,npg)
  call alloc_array2d(bvxd_bh,0,ntg,0,npg)
  call alloc_array2d(bvyd_bh,0,ntg,0,npg)
  call alloc_array2d(bvzd_bh,0,ntg,0,npg)
  call alloc_array2d(sous1, 0,ntg,0,npg)
  call alloc_array2d(sous2, 0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf ,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
!
  psi_bh(0:ntg,0:npg)  = psi(0,0:ntg,0:npg)
  bvxd_bh(0:ntg,0:npg) = bvxd(0,0:ntg,0:npg)
  bvyd_bh(0:ntg,0:npg) = bvyd(0,0:ntg,0:npg)
  bvzd_bh(0:ntg,0:npg) = bvzd(0,0:ntg,0:npg)
!
  call sourceterm_surface_int(alph,0,sou_bhsurf,dsou_bhsurf)
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(val,psi_bh,itg,ipg)
      psi2 = val**2
      psi6 = val**6
      sous1(itg,ipg) = -psi2*dsou_bhsurf(itg,ipg)
!
      call interpo_linear_type0_2Dsurf(bv(1),bvxd_bh,itg,ipg)
      call interpo_linear_type0_2Dsurf(bv(2),bvyd_bh,itg,ipg)
      call interpo_linear_type0_2Dsurf(bv(3),bvzd_bh,itg,ipg)
      sous2(itg,ipg)=0.0d0
      do ib = 1, 3
        do ia = 1, 3
          work(1:2,1:2) = tfkij_grid(0,itg-1:itg,ipg-1:ipg,ia,ib)
          call interpo_linear1p_type0_2Dsurf(Aij_surf,work)
!         Aij_surf = tfkij(1,itg,ipg,ia,ib)
!
          ni = -hvec_bh_cbh_xg(itg,ipg,ib)/rg(0)
          vphi_cm = hvec_bh_cm_phig(itg,ipg,ia)
          beta = bv(ia) + ome_bh*vphi_cm
!
          sous2(itg,ipg) = sous2(itg,ipg) + Aij_surf*beta*ni
        end do
      end do
      sous2(itg,ipg) = sous2(itg,ipg)*psi6
!
      sous(itg,ipg)  = sous1(itg,ipg) - sous2(itg,ipg)
    end do
  end do
!
  deallocate(psi_bh)
  deallocate(bvxd_bh)
  deallocate(bvyd_bh)
  deallocate(bvzd_bh)
  deallocate(sous1)
  deallocate(sous2)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
end subroutine source_ang_mom_smarr

