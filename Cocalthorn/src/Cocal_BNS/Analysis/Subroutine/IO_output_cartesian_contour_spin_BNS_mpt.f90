subroutine IO_output_cartesian_contour_spin_BNS_mpt(impt)
  use phys_constant, only : long
  use def_metric_cartesian
  use def_matter_cartesian
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: emdc, vxc, vyc, vzc, psic, alphc, bvxc, bvyc, bvzc, vepc
  real(long) :: hhc, prec, rhoc, enec
  real(long) :: wxc, wyc, wzc
  integer :: ix, iy, iz, impt
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/), char_1
!
  open(14,file='bns_contour_xy_mpt'//np(impt)//'.dat',status='unknown')
  iz = nx_mid
  do iy = 1, ny
    write(14,*) ' '
    do ix = 1, nx
      emdc = emdca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      psic = psica(ix,iy,iz)
      alphc=alphca(ix,iy,iz)
      bvxc = bvxdca(ix,iy,iz)
      bvyc = bvydca(ix,iy,iz)
      bvzc = bvzdca(ix,iy,iz)
      vepc = vepca(ix,iy,iz)
      wxc  = wxspca(ix,iy,iz)
      wyc  = wyspca(ix,iy,iz)
      wzc  = wzspca(ix,iy,iz)
      if (dabs(emdc).le.small) emdc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxc).le.small) bvxc = 0.0d0
      if (dabs(bvyc).le.small) bvyc = 0.0d0
      if (dabs(bvzc).le.small) bvzc = 0.0d0
      if (dabs(vepc).le.small) vepc = 0.0d0
      if (dabs(wxc).le.small) wxc = 0.0d0
      if (dabs(wyc).le.small) wyc = 0.0d0
      if (dabs(wzc).le.small) wzc = 0.0d0
!
      call peos_q2hprho(emdc, hhc, prec, rhoc, enec)

      write(14,'(20es14.6)') x(ix), y(iy), emdc, vxc, vyc, vzc,  &
     &  psic, alphc,bvxc,bvyc,bvzc,vepc,wxc,wyc,wzc, hhc, rhoc
!
    end do
  end do
  close(14)
!
  open(14,file='bns_contour_xz_mpt'//np(impt)//'.dat',status='unknown')
  iy = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do ix = 1, nx
      emdc = emdca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      psic = psica(ix,iy,iz)
      alphc=alphca(ix,iy,iz)
      bvxc = bvxdca(ix,iy,iz)
      bvyc = bvydca(ix,iy,iz)
      bvzc = bvzdca(ix,iy,iz)
      vepc = vepca(ix,iy,iz)
      wxc  = wxspca(ix,iy,iz)
      wyc  = wyspca(ix,iy,iz)
      wzc  = wzspca(ix,iy,iz)
      if (dabs(emdc).le.small) emdc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxc).le.small) bvxc = 0.0d0
      if (dabs(bvyc).le.small) bvyc = 0.0d0
      if (dabs(bvzc).le.small) bvzc = 0.0d0
      if (dabs(vepc).le.small) vepc = 0.0d0
      if (dabs(wxc).le.small) wxc = 0.0d0
      if (dabs(wyc).le.small) wyc = 0.0d0
      if (dabs(wzc).le.small) wzc = 0.0d0
!
      call peos_q2hprho(emdc, hhc, prec, rhoc, enec)

      write(14,'(20es14.6)') x(ix), z(iz),  emdc, vxc, vyc, vzc,  &
     &  psic, alphc, bvxc, bvyc, bvzc, vepc,wxc,wyc,wzc, hhc, rhoc
!
    end do
  end do
  close(14)
!
  open(14,file='bns_contour_yz_mpt'//np(impt)//'.dat',status='unknown')
  ix = nx_mid
  do iz = 1, nz
    write(14,*) ' '
    do iy = 1, ny
      emdc = emdca(ix,iy,iz)
      vxc  =  vxca(ix,iy,iz)
      vyc  =  vyca(ix,iy,iz)
      vzc  =  vzca(ix,iy,iz)
      psic = psica(ix,iy,iz)
      alphc=alphca(ix,iy,iz)
      bvxc = bvxdca(ix,iy,iz)
      bvyc = bvydca(ix,iy,iz)
      bvzc = bvzdca(ix,iy,iz)
      vepc = vepca(ix,iy,iz)
      wxc  = wxspca(ix,iy,iz)
      wyc  = wyspca(ix,iy,iz)
      wzc  = wzspca(ix,iy,iz)
      if (dabs(emdc).le.small) emdc = 0.0d0
      if (dabs(vxc).le.small) vxc = 0.0d0
      if (dabs(vyc).le.small) vyc = 0.0d0
      if (dabs(vzc).le.small) vzc = 0.0d0
      if (dabs(psic).le.small) psic = 0.0d0
      if (dabs(alphc).le.small) alphc = 0.0d0
      if (dabs(bvxc).le.small) bvxc = 0.0d0
      if (dabs(bvyc).le.small) bvyc = 0.0d0
      if (dabs(bvzc).le.small) bvzc = 0.0d0
      if (dabs(vepc).le.small) vepc = 0.0d0
      if (dabs(wxc).le.small) wxc = 0.0d0
      if (dabs(wyc).le.small) wyc = 0.0d0
      if (dabs(wzc).le.small) wzc = 0.0d0
!
      call peos_q2hprho(emdc, hhc, prec, rhoc, enec)

      write(14,'(20es14.6)') y(iy), z(iz), emdc, vxc, vyc, vzc,  &
     &  psic, alphc, bvxc, bvyc, bvzc, vepc,wxc,wyc,wzc, hhc, rhoc
!
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_spin_BNS_mpt
