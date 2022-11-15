subroutine bh_boundary_AH_n_alps(dsou_surf)
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
!
  do ipg = 1, npg
    do itg = 1, ntg
      dsou_surf(itg,ipg)= 0.5d0
    end do
  end do
!
end subroutine bh_boundary_AH_n_alps
