subroutine IO_output_cartesian_contour_potential_test_mpt(impt)
  use phys_constant, only : long
  use def_metric_cartesian, only  : psica
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: psic
  integer :: ix, iy, iz
  integer :: impt
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/)
!
  open(14,file='rns_contour_xy_mpt'//np(impt)//'.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
      psic = psica(ix,iy,iz)
      if (dabs(psic).le.small) psic = 0.0d0
      write(14,'(20es14.6)') x(ix), y(iy), &
     &  psic
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_xz_mpt'//np(impt)//'.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
      psic = psica(ix,iy,iz)
      if (dabs(psic).le.small) psic = 0.0d0
      write(14,'(20es14.6)') x(ix), z(iz), &
     &  psic
!
    end do
  end do
  close(14)
!
  open(14,file='rns_contour_yz_mpt'//np(impt)//'.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
      psic = psica(ix,iy,iz)
      if (dabs(psic).le.small) psic = 0.0d0
      write(14,'(20es14.6)') y(iy), z(iz), &
     &  psic
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_potential_test_mpt
