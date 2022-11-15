subroutine IO_output_cartesian_contour_CF
  use phys_constant, only : long
  use def_metric_cartesian
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: psic, alphc, bvxdc, bvydc, bvzdc
  integer :: ix, iy, iz
!
  write(6,*) "Output rns_contour_**_CF.dat..."

  open(14,file='rns_contour_xy_CF.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
      psic = psica(ix,iy,iz)
      alphc = alphca(ix,iy,iz)
      bvxdc = bvxdca(ix,iy,iz)
      bvydc = bvydca(ix,iy,iz)
      bvzdc = bvzdca(ix,iy,iz)
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxdc).le.small) bvxdc = 0.0d0
      if (dabs(bvydc).le.small) bvydc = 0.0d0
      if (dabs(bvzdc).le.small) bvzdc = 0.0d0
      write(14,'(20es14.6)') x(ix), y(iy), psic, alphc, bvxdc, bvydc, bvzdc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_xz_CF.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
      psic = psica(ix,iy,iz)
      alphc = alphca(ix,iy,iz)
      bvxdc = bvxdca(ix,iy,iz)
      bvydc = bvydca(ix,iy,iz)
      bvzdc = bvzdca(ix,iy,iz)
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxdc).le.small) bvxdc = 0.0d0
      if (dabs(bvydc).le.small) bvydc = 0.0d0
      if (dabs(bvzdc).le.small) bvzdc = 0.0d0
      write(14,'(20es14.6)') x(ix), z(iz), psic, alphc, bvxdc, bvydc, bvzdc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_yz_CF.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
      psic = psica(ix,iy,iz)
      alphc = alphca(ix,iy,iz)
      bvxdc = bvxdca(ix,iy,iz)
      bvydc = bvydca(ix,iy,iz)
      bvzdc = bvzdca(ix,iy,iz)
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxdc).le.small) bvxdc = 0.0d0
      if (dabs(bvydc).le.small) bvydc = 0.0d0
      if (dabs(bvzdc).le.small) bvzdc = 0.0d0
      write(14,'(20es14.6)') y(iy), z(iz), psic, alphc, bvxdc, bvydc, bvzdc
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_CF
