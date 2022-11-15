subroutine copy_hgfn_nb_to_hgfn
  use grid_parameter, only : nrg, nlg
  use copy_array_3d
  use radial_green_fn_grav, only : hgfn, gfnsf
  use radial_green_fn_grav_bhex_nb, only : hgfn_nb, gfnsf_nb
  implicit none
  call copy_array3d(hgfn_nb,hgfn,1,nrg,0,nlg,0,nrg)
  call copy_array3d(gfnsf_nb,gfnsf,0,nlg,0,nrg,1,4)
end subroutine copy_hgfn_nb_to_hgfn
