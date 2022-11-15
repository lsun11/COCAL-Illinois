subroutine copy_to_hgfn_and_gfnsf(hgfn_bc,gfnsf_bc)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, nlg
  use copy_array_3d
  use radial_green_fn_grav, only : hgfn, gfnsf
  implicit none
  real(long), pointer :: hgfn_bc(:,:,:), gfnsf_bc(:,:,:)
  call copy_array3d(hgfn_bc,hgfn,1,nrg,0,nlg,0,nrg)
  call copy_array3d(gfnsf_bc,gfnsf,0,nlg,0,nrg,1,4)
end subroutine copy_to_hgfn_and_gfnsf
