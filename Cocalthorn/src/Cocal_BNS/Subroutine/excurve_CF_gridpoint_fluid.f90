subroutine excurve_CF_gridpoint_fluid
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_metric_excurve_grid, only : tfkij_grid
  use def_excurve_on_SFC_CF
  use interface_interpo_gr2fl
  use make_array_3d
  implicit none
  integer :: ia, ib
  real(long), pointer :: tempf(:,:,:), tempg(:,:,:)

  call alloc_array3d(tempf,0,nrf,0,ntf,0,npf)
  call alloc_array3d(tempg,0,nrg,0,ntg,0,npg)
!
  do ib = 1, 3
    do ia = 1, 3
      tempg(0:nrg,0:ntg,0:npg) = tfkij_grid(0:nrg,0:ntg,0:npg,ia,ib)
      call interpo_gr2fl(tempg, tempf)
      tfkij_grid_fluid(0:nrf,0:ntf,0:npf,ia,ib) = tempf(0:nrf,0:ntf,0:npf)
    end do
  end do
!
  deallocate(tempf)
  deallocate(tempg)
  
end subroutine excurve_CF_gridpoint_fluid
