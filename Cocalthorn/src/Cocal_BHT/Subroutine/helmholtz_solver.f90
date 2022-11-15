subroutine helmholtz_solver(sou,pot)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, npg, ntg
  use def_matter_parameter, only : ome
  use radial_green_fn_hrethadv
  use make_array_3d
  use interface_copy_to_bsjy_and_sbsjy
  use interface_helmholtz_solver_vol_int
  implicit none
  real(long), pointer :: pot(:,:,:), sou(:,:,:)
  real(long), pointer :: pot_vol(:,:,:)
! 
  call alloc_array3d(pot_vol,0,nrg,0,ntg,0,npg)
  call calc_radial_green_fn_hrethadv(ome)
  call copy_to_bsjy_and_sbsjy(bsjy_hrha,sbsjy_hrha,sbsjyp_hrha)
!
  call helmholtz_solver_vol_int(sou,pot_vol)
  pot(0:nrg,0:ntg,0:npg) = &
  & pot_vol(0:nrg,0:ntg,0:npg)
!
  deallocate(pot_vol)
! 
end subroutine helmholtz_solver
