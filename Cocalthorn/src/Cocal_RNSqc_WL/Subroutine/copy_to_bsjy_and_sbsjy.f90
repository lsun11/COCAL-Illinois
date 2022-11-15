subroutine copy_to_bsjy_and_sbsjy(bsjy_bc,sbsjy_bc,sbsjyp_bc)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, nlg
  use copy_array_4d
  use copy_array_3d
  use radial_green_fn_helmholtz
  implicit none
  real(long), pointer :: bsjy_bc(:,:,:,:)
  real(long), pointer :: sbsjy_bc(:,:,:), sbsjyp_bc(:,:,:)
  call copy_array4d(bsjy_bc,bsjy,1,nrg,0,nlg,0,nlg,0,nrg)
  call copy_array3d(sbsjy_bc,sbsjy,0,nlg,0,nlg,0,nrg)
  call copy_array3d(sbsjyp_bc,sbsjyp,0,nlg,0,nlg,0,nrg)
end subroutine copy_to_bsjy_and_sbsjy
