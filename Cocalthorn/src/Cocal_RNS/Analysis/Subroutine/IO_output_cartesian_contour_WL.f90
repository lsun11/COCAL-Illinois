subroutine IO_output_cartesian_contour_WL
  use phys_constant, only : long
  use def_metric_hij_cartesian, only  : hxxdca, hxydca, hxzdca, &
  &                                     hyydca, hyzdca, hzzdca
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: hxxdc, hxydc, hxzdc, hyydc, hyzdc, hzzdc
  integer :: ix, iy, iz
!
  write(6,*)  "Output rns_contour_hij_**.dat..."

  open(14,file='rns_contour_hij_xy.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
!
      hxxdc = hxxdca(ix,iy,iz)
      hxydc = hxydca(ix,iy,iz)
      hxzdc = hxzdca(ix,iy,iz)
      hyydc = hyydca(ix,iy,iz)
      hyzdc = hyzdca(ix,iy,iz)
      hzzdc = hzzdca(ix,iy,iz)
      if (dabs(hxxdc).le.small) hxxdc = 0.0d0
      if (dabs(hxydc).le.small) hxydc = 0.0d0
      if (dabs(hxzdc).le.small) hxzdc = 0.0d0
      if (dabs(hyydc).le.small) hyydc = 0.0d0
      if (dabs(hyzdc).le.small) hyzdc = 0.0d0
      if (dabs(hzzdc).le.small) hzzdc = 0.0d0
!
      write(14,'(20es14.6)') x(ix), y(iy), &
     &  hxxdc, hxydc, hxzdc, hyydc, hyzdc, hzzdc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_hij_xz.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
!
      hxxdc = hxxdca(ix,iy,iz)
      hxydc = hxydca(ix,iy,iz)
      hxzdc = hxzdca(ix,iy,iz)
      hyydc = hyydca(ix,iy,iz)
      hyzdc = hyzdca(ix,iy,iz)
      hzzdc = hzzdca(ix,iy,iz)
      if (dabs(hxxdc).le.small) hxxdc = 0.0d0
      if (dabs(hxydc).le.small) hxydc = 0.0d0
      if (dabs(hxzdc).le.small) hxzdc = 0.0d0
      if (dabs(hyydc).le.small) hyydc = 0.0d0
      if (dabs(hyzdc).le.small) hyzdc = 0.0d0
      if (dabs(hzzdc).le.small) hzzdc = 0.0d0
!
      write(14,'(20es14.6)') x(ix), z(iz), &
     &  hxxdc, hxydc, hxzdc, hyydc, hyzdc, hzzdc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_hij_yz.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
!
      hxxdc = hxxdca(ix,iy,iz)
      hxydc = hxydca(ix,iy,iz)
      hxzdc = hxzdca(ix,iy,iz)
      hyydc = hyydca(ix,iy,iz)
      hyzdc = hyzdca(ix,iy,iz)
      hzzdc = hzzdca(ix,iy,iz)
      if (dabs(hxxdc).le.small) hxxdc = 0.0d0
      if (dabs(hxydc).le.small) hxydc = 0.0d0
      if (dabs(hxzdc).le.small) hxzdc = 0.0d0
      if (dabs(hyydc).le.small) hyydc = 0.0d0
      if (dabs(hyzdc).le.small) hyzdc = 0.0d0
      if (dabs(hzzdc).le.small) hzzdc = 0.0d0
!
      write(14,'(20es14.6)') y(iy), z(iz), &
     &  hxxdc, hxydc, hxzdc, hyydc, hyzdc, hzzdc
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_WL
