subroutine copy_hgfn_di_to_hgfn
  use grid_parameter, only : nrg, nlg
  use copy_array_3d
  use radial_green_fn_grav, only : hgfn, gfnsf
  use radial_green_fn_grav_bhex_di, only : hgfn_di, gfnsf_di
  implicit none
  call copy_array3d(hgfn_di,hgfn,1,nrg,0,nlg,0,nrg)
  call copy_array3d(gfnsf_di,gfnsf,0,nlg,0,nrg,1,4)
end subroutine copy_hgfn_di_to_hgfn
