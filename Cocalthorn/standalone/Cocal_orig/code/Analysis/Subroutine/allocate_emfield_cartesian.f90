subroutine allocate_emfield_cartesian
  use phys_constant, only : long
  use grid_parameter_cartesian, only : nx, ny, nz
  use def_emfield_cartesian
  use def_faraday_tensor_cartesian
  use make_array_3d
  use make_array_4d
  implicit none
!
  call alloc_array3d(vaca,1,nx,1,ny,1,nz)
  call alloc_array3d(vaxdca,1,nx,1,ny,1,nz)
  call alloc_array3d(vaydca,1,nx,1,ny,1,nz)
  call alloc_array3d(vazdca,1,nx,1,ny,1,nz)
  call alloc_array3d(fxd_gridca,1,nx,1,ny,1,nz)
  call alloc_array3d(fyd_gridca,1,nx,1,ny,1,nz)
  call alloc_array3d(fzd_gridca,1,nx,1,ny,1,nz)
  call alloc_array4d(fijd_gridca,1,nx,1,ny,1,nz,1,3)
!
end subroutine allocate_emfield_cartesian
