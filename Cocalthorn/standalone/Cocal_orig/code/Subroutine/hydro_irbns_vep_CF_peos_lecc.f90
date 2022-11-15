subroutine hydro_irbns_vep_CF_peos_lecc(vpot,iter,impt,ihy)
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use coordinate_grav_r, only : rg, hrg
  use def_matter, only : rs
  use make_array_2d
  use make_array_3d
  use weight_midpoint_fluid_sphcoord
  use interface_source_vep_CF_peos_lecc
!  use interface_source_vep_CF_peos_v1
  use interface_interpo_flsfc2flsph_midpoint
  use interface_poisson_solver_fluid_sphcoord
  use interface_interpo_flsph2flsfc
  use interface_source_vep_surface_CF_peos
  use interface_poisson_solver_homogeneous_sol_lecc
  use radial_green_fn_grav

  use interface_interpo_linear_type0_2Dsurf

!  use interface_copy_to_hgfn_and_gfnsf

  implicit none
  real(long), pointer :: vpot(:,:,:)
  real(long), pointer :: sou(:,:,:), soufc(:,:,:), surp(:,:)
  real(long), pointer :: vpotfc(:,:,:), vpot_v(:,:,:), vpot_b(:,:,:)
  integer :: ir,it,ip, iter, impt, ihy
  real(long) :: rrff
  character(30) :: char1, char2, char3, char4, char5, char6, char7
!
  call alloc_array3d(sou,    0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(soufc,  0, nrf, 0, ntf, 0, npf)
  call alloc_array2d(surp,   0, ntf, 0, npf)
  call alloc_array3d(vpotfc, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vpot_v, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(vpot_b, 0, nrf, 0, ntf, 0, npf)
!
  call source_vep_CF_peos_lecc(sou)

  call interpo_flsfc2flsph_midpoint(sou,soufc)

  call calc_weight_midpoint_fluid_sphcoord

  call calc_hgfn

  call poisson_solver_fluid_sphcoord(soufc,vpotfc)

  call interpo_flsph2flsfc(vpotfc,vpot_v)

  call source_vep_surface_CF_peos(vpot_v,surp)

  call poisson_solver_homogeneous_sol_lecc(surp,vpot_b)

  if (mod(iter,5)==0.and.impt==1.and.ihy==4) then
    write(char1, '(i5)') iter
    char2 = adjustl(char1)
    write(char4, '(i5)') impt
    char5 = adjustl(char4)
    write(char6, '(i5)') ihy
    char7 = adjustl(char6)

!    char3 = 'iteration' // trim(char2) // '_mpt' // trim(char5) // '_ihy' // trim(char7) //'.txt'
    char3 = 'iteration' // trim(char2) // '_mpt' // trim(char5) // '.txt'
    open(12,file=char3,status='unknown')
    it=12;  ip=6
    call interpo_linear_type0_2Dsurf(rrff,rs,it,ip)
    write(12,'(a1,a2,4a16,a30)') '#', 'ir', 'hrg*rrff', 'sou', 'hrg', 'soufc', 'plot using every :::0::0'
    do ir=1,nrf
      write(12,'(i3,4e16.6)')  ir, hrg(ir)*rrff, sou(ir,it,ip),  hrg(ir), soufc(ir,it,ip)
    end do
    write(12,*) '#------------------------------------------------------------------------------------------------------'
    write(12,*) ""
    write(12,'(a1,a2,5a16,a30)') '#', 'ir', 'rg', 'vpotfc', 'rs*rg', 'vpot_v', 'vpot_b', 'plot using every :::1::1'
    do ir=0,nrf
      write(12,'(i3,5e16.6)')  ir, rg(ir), vpotfc(ir,it,ip),  &
        &            rs(it,ip)*rg(ir), vpot_v(ir,it,ip), vpot_b(ir,it,ip)
    end do
    write(12,*) ""
    close(12)
  endif

  vpot(0:nrf,0:ntf,0:npf) = vpot_v(0:nrf,0:ntf,0:npf) + vpot_b(0:nrf,0:ntf,0:npf)
!
  deallocate(sou)
  deallocate(soufc)
  deallocate(surp)
  deallocate(vpotfc)
  deallocate(vpot_v)
  deallocate(vpot_b)
end subroutine hydro_irbns_vep_CF_peos_lecc
