subroutine hydro_irbns_vep_WL_peos(vpot)
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use make_array_2d
  use make_array_3d
  implicit none
  real(long), pointer :: vpot(:,:,:)
  real(long), pointer :: sou(:,:,:), soufc(:,:,:)
  real(long), pointer :: vpotfc(:,:,:), vpot_v(:,:,:), vpot_b(:,:,:)
!
  call alloc_array3d(sou,    0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(soufc,  0, nrf, 0, ntf, 0, npf)
  call alloc_array2d(surp,   0, ntf, 0, npf)
  call alloc_array3d(vpotfc, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vpot_v, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vpot_b, 0, nrf, 0, ntf, 0, npf)
!
  call source_vep_WL_peos(sou)
  call interpo_flsfc2flsph_midpoint(sou,soufc)
!
  call calc_weight_midpoint_fluid_sphcoord
  call poisson_solver_fluid_sphcoord(soufc,vpotfc)
  call interpo_flsph2flsfc(vpotfc,vpot_v)
!
  call source_vep_surface_WL_peos(vpot_v,surp)
  call poisson_solver_homogeneous_sol(surp,vpot_b)
  vpot(0:nrf,0:ntf,0:npf) = vpot_v(0:nrf,0:ntf,0:npf)
  &                       + vpot_b(0:nrf,0:ntf,0:npf)
!
  deallocate(sou)
  deallocate(soufc)
  deallocate(surp)
  deallocate(vpotfc)
  deallocate(vpot_v)
  deallocate(vpot_b)
end subroutine hydro_irbns_vep_WL_peos
