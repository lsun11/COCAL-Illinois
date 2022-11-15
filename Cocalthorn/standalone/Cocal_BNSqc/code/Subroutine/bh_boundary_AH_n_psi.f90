subroutine bh_boundary_AH_n_psi(dsou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use def_metric, only : psi
  use def_bh_parameter, only : ome_bh, spin_bh, alph_bh, psi_bh, bh_bctype
  use def_binary_parameter,    only : sepa
  use def_metric_excurve_grid, only : tfkij_grid
  use def_vector_bh, only : hvec_bh_cbh_xg,hvec_bh_cm_phig,hvec_bh_cbh_spin
  use make_array_2d
  use interface_grdr_gridpoint_type0_nosym
  use interface_interpo_linear_type0_2Dsurf
  implicit none
  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long), pointer :: psi_bhsurf(:,:), dpsi_bhsurf(:,:)
  real(long) :: work(2,2), deriv, dpsidr
  real(long) :: Aij_surf, Aij_sisj, si, sj
  real(long) :: psigc, psigc3, psigcinv2
  real(long) :: vphi_cm, vphi_cbh
  integer    :: itg, ipg, ii, jj

  call alloc_array2d(psi_bhsurf, 0, ntg, 0, npg)
!
  psi_bhsurf(0:ntg,0:npg) = psi(0,0:ntg,0:npg)
!
  do ipg = 1, npg
    do itg = 1, ntg
      call interpo_linear_type0_2Dsurf(psigc,psi_bhsurf,itg,ipg)
      psigc3 = psigc**3
      Aij_sisj = 0.0d0
      do ii = 1, 3
        do jj = 1, 3
          work(1:2,1:2) = tfkij_grid(0,itg-1:itg,ipg-1:ipg,ii,jj)
          call interpo_linear1p_type0_2Dsurf(Aij_surf,work)
          si = hvec_bh_cbh_xg(itg,ipg,ii)/rgin
          sj = hvec_bh_cbh_xg(itg,ipg,jj)/rgin
          Aij_sisj = Aij_sisj + Aij_surf*si*sj
        end do
      end do
      dsou_surf(itg,ipg)= - 0.5d0*psigc/rgin  - 0.25d0*psigc3*Aij_sisj
    end do
  end do

  deallocate(psi_bhsurf)
!
end subroutine bh_boundary_AH_n_psi
