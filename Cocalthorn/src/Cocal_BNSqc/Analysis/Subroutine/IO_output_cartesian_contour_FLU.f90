subroutine IO_output_cartesian_contour_FLU
  use phys_constant, only : long
  use def_matter_cartesian, only  : emdca, vxca, vyca, vzca, utca
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: emdc, vxc, vyc, vzc, utc
  integer :: ix, iy, iz
!
  open(14,file='rns_contour_xy_FLU.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
      emdc = emdca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      utc  =  utca(ix,iy,iz)
      if (dabs(emdc).le.small) emdc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(utc).le.small) utc = 0.0d0
!
      write(14,'(20es14.6)') x(ix), y(iy), emdc, vxc, vyc, vzc, utc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_xz_FLU.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
      emdc = emdca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      utc  =  utca(ix,iy,iz)
      if (dabs(emdc).le.small) emdc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(utc).le.small) utc = 0.0d0
!
      write(14,'(20es14.6)') x(ix), z(iz), emdc, vxc, vyc, vzc, utc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_yz_FLU.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
      emdc = emdca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      utc  =  utca(ix,iy,iz)
      if (dabs(emdc).le.small) emdc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(utc).le.small) utc = 0.0d0
!
      write(14,'(20es14.6)') y(iy), z(iz), emdc, vxc, vyc, vzc, utc
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_FLU
