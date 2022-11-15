subroutine reset_bh_boundary_AH
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use def_metric
  use def_bh_parameter, only : ome_bh, spin_bh, alph_bh
  use def_binary_parameter,    only : sepa
  use def_metric_excurve_grid, only : tfkij_grid
  use def_vector_bh, only : vec_bh_cbh_xg, vec_bh_cm_phig, vec_bh_cbh_phig
  use make_array_2d
  use make_array_3d
  use interface_grdr_gridpoint_type0_nosym
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long), pointer :: psi_bhsurf(:,:), dpsi_bhsurf(:,:)
  real(long), pointer :: bvxyzd(:,:,:)
  real(long) :: work(2,2), deriv, dpsidr
  real(long) :: Aij_surf, Aij_sisj, si, sj
  real(long) :: psigc, psigc3, psigcinv2
  real(long) :: vphi_cm, vphi_cbh
  integer    :: itg, ipg, ii, jj
!
  call alloc_array2d(psi_bhsurf, 0, ntg, 0, npg)
  call alloc_array2d(dpsi_bhsurf, 0, ntg, 0, npg)
  call alloc_array3d(bvxyzd,0,ntg,0,npg,1,3)
!
  psi_bhsurf(0:ntg,0:npg) = psi(0,0:ntg,0:npg)
!
  do ipg = 0, npg
    do itg = 0, ntg
      call grdr_gridpoint_type0_nosym(psi,deriv,0,itg,ipg)
      dpsi_bhsurf(itg,ipg) = deriv
    end do
  end do
!
  do ipg = 0, npg
    do itg = 0, ntg
      psigc = psi(0,itg,ipg)
      psigc3 = psigc**3
      Aij_sisj = 0.0d0
      do ii = 1, 3
        do jj = 1, 3
          Aij_surf = tfkij_grid(0,itg,ipg,ii,jj)
          si = vec_bh_cbh_xg(itg,ipg,ii)/rgin
          sj = vec_bh_cbh_xg(itg,ipg,jj)/rgin
          Aij_sisj = Aij_sisj + Aij_surf*si*sj
        end do
      end do
      psi(0,itg,ipg) = 2.0d0*rgin*(-dpsi_bhsurf(itg,ipg) - 0.25d0*psigc3*Aij_sisj)
!      dsou_surf(itg,ipg)= - 0.5d0*psigc/rgin  - 0.25d0*psigc3*Aij_sisj
    end do
  end do
  alps(0,0:ntg,0:npg) = alph_bh*psi(0,0:ntg,0:npg)
!  dsou_surf(1:ntg,1:npg)= alph_bh*dsou_surf(1:ntg,1:npg)
!
  do ii = 1,3
    do ipg = 0, npg
      do itg = 0, ntg
        psigc = psi(0,itg,ipg)
        psigcinv2 = 1.0d0/psigc**2
        si = vec_bh_cbh_xg(itg,ipg,ii)/rgin
        vphi_cm  = vec_bh_cm_phig(itg,ipg,ii)
        vphi_cbh = vec_bh_cbh_phig(itg,ipg,ii)
        bvxyzd(itg,ipg,ii) = alph_bh*psigcinv2*si - ome_bh*vphi_cm - spin_bh*vphi_cbh
      end do
    end do
  end do
!
  bvxd(0,0:ntg,0:npg) = bvxyzd(0:ntg,0:npg,1)
  bvyd(0,0:ntg,0:npg) = bvxyzd(0:ntg,0:npg,2)
  bvzd(0,0:ntg,0:npg) = bvxyzd(0:ntg,0:npg,3)
!
  deallocate(bvxyzd)
  deallocate(psi_bhsurf)
  deallocate(dpsi_bhsurf)
end subroutine reset_bh_boundary_AH
