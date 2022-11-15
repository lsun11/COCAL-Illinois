subroutine IO_output_cartesian_contour_qeos
  use phys_constant, only : long
  use def_metric_cartesian, only  : alphca, psica, bvxdca, bvydca, bvzdca
  use def_matter_parameter, only  : rhos_qs
  use def_matter_cartesian_qeos, only  : rhoca, vxca, vyca, vzca
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z

  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: rhoc, vxc, vyc, vzc, psic, alphc, bvxc, bvyc, bvzc
  integer :: ix, iy, iz
!
  open(14,file='rns_contour_xy.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
      rhoc = rhoca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      psic = psica(ix,iy,iz)
      alphc=alphca(ix,iy,iz)
      bvxc = bvxdca(ix,iy,iz)
      bvyc = bvydca(ix,iy,iz)
      bvzc = bvzdca(ix,iy,iz)
      if (dabs(rhoc).le.rhos_qs) rhoc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxc).le.small) bvxc = 0.0d0
      if (dabs(bvyc).le.small) bvyc = 0.0d0
      if (dabs(bvzc).le.small) bvzc = 0.0d0
!
      write(14,'(20es14.6)') x(ix), y(iy), &
     &  rhoc, vxc, vyc, vzc, psic, alphc, bvxc, bvyc, bvzc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_xz.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
      rhoc = rhoca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      psic = psica(ix,iy,iz)
      alphc=alphca(ix,iy,iz)
      bvxc = bvxdca(ix,iy,iz)
      bvyc = bvydca(ix,iy,iz)
      bvzc = bvzdca(ix,iy,iz)
      if (dabs(rhoc).le.rhos_qs) rhoc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxc).le.small) bvxc = 0.0d0
      if (dabs(bvyc).le.small) bvyc = 0.0d0
      if (dabs(bvzc).le.small) bvzc = 0.0d0
!
      write(14,'(20es14.6)') x(ix), z(iz), &
     &  rhoc, vxc, vyc, vzc, psic, alphc, bvxc, bvyc, bvzc
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_yz.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
      rhoc = rhoca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      psic = psica(ix,iy,iz)
      alphc=alphca(ix,iy,iz)
      bvxc = bvxdca(ix,iy,iz)
      bvyc = bvydca(ix,iy,iz)
      bvzc = bvzdca(ix,iy,iz)
      if (dabs(rhoc).le.rhos_qs) rhoc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxc).le.small) bvxc = 0.0d0
      if (dabs(bvyc).le.small) bvyc = 0.0d0
      if (dabs(bvzc).le.small) bvzc = 0.0d0
!
      write(14,'(20es14.6)') y(iy), z(iz), &
     &  rhoc, vxc, vyc, vzc, psic, alphc, bvxc, bvyc, bvzc
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_qeos
