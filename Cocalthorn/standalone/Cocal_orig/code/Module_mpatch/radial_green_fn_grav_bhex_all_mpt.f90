! --- Computation of the radial green's function hgfn. ---
! ________________________________________________________
module radial_green_fn_grav_bhex_all_mpt
  use phys_constant, only : long
  use grid_parameter, only : nrg, nlg, rgin, rgout
  use coordinate_grav_r, only : rg, rginv, hrg, hrginv
  use make_array_3d
  implicit none
  real(long), pointer  ::  hgfn_nb_(:,:,:,:), gfnsf_nb_(:,:,:,:)
  real(long), pointer  ::  hgfn_di_(:,:,:,:), gfnsf_di_(:,:,:,:)
  real(long), pointer  ::  hgfn_dd_(:,:,:,:), gfnsf_dd_(:,:,:,:)
  real(long), pointer  ::  hgfn_nd_(:,:,:,:), gfnsf_nd_(:,:,:,:)
end module radial_green_fn_grav_bhex_all_mpt
