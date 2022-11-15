subroutine allocate_metric_and_matter_WL_cartesian
  use phys_constant, only : long
  use grid_parameter_cartesian, only : nx, ny, nz
  use def_metric_cartesian
  use def_metric_hij_cartesian
  use def_matter_cartesian
  use make_array_3d
  implicit none
!
  call alloc_array3d(emdca,1,nx,1,ny,1,nz)
  call alloc_array3d(omeca,1,nx,1,ny,1,nz)
  call alloc_array3d(utca,1,nx,1,ny,1,nz)
  call alloc_array3d(vxca,1,nx,1,ny,1,nz)
  call alloc_array3d(vyca,1,nx,1,ny,1,nz)
  call alloc_array3d(vzca,1,nx,1,ny,1,nz)
  call alloc_array3d(psica,1,nx,1,ny,1,nz)
  call alloc_array3d(alphca,1,nx,1,ny,1,nz)
  call alloc_array3d(bvxdca,1,nx,1,ny,1,nz)
  call alloc_array3d(bvydca,1,nx,1,ny,1,nz)
  call alloc_array3d(bvzdca,1,nx,1,ny,1,nz)
  call alloc_array3d(bvxuca,1,nx,1,ny,1,nz)
  call alloc_array3d(bvyuca,1,nx,1,ny,1,nz)
  call alloc_array3d(bvzuca,1,nx,1,ny,1,nz)
  call alloc_array3d(hxxdca,1,nx,1,ny,1,nz)
  call alloc_array3d(hxydca,1,nx,1,ny,1,nz)
  call alloc_array3d(hxzdca,1,nx,1,ny,1,nz)
  call alloc_array3d(hyydca,1,nx,1,ny,1,nz)
  call alloc_array3d(hyzdca,1,nx,1,ny,1,nz)
  call alloc_array3d(hzzdca,1,nx,1,ny,1,nz)
!
end subroutine allocate_metric_and_matter_WL_cartesian
