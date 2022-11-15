!  Cartesian_coordinate
!______________________________________________
module coordinate_grav_xyz
  use grid_parameter_cartesian, only : nx, ny, nz, nxhalf, nnx, nny, nnz, nstar
  implicit none
  real(8) :: dx, dy, dz
  real(8) :: x(nnx), y(nny), z(nnz)
! Subroutine
contains
!
subroutine grid_xyz
  implicit none
  Integer  ::  ix, iy, iz, igrid
!
  igrid = mod(nx,2)
  if (igrid.eq.0) write(6,*) ' ## Volume centered ## '
  if (igrid.eq.1) write(6,*) ' ## Grid centered ## '
!
  dx = 1.0d0/dble(nstar)
  dy = dx
  dz = dx
! Origines of Cartesian coordinate and the spherical coordinate agree.
  x(1) = - dx*dble(nxhalf) + 0.5d0*dx*dble(1-igrid)
  y(1) = x(1)
  z(1) = x(1)
  do ix = 2, nx
    x(ix) = x(1) + dx*dble(ix-1)
    y(ix) = x(ix)
    z(ix) = x(ix)
  end do
! 
end subroutine grid_xyz
end module coordinate_grav_xyz
