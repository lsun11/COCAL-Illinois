subroutine IO_output_cartesian_planes
  use phys_constant, only : long
  use def_metric, only : bvxd, bvyd, bvzd
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use coordinate_grav_xyz, only  : x, y, z
  use def_binary_parameter,    only : sepa, dis
  use grid_parameter, only : rgin
  use interface_modules_cartesian
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: psic, alphc, bvxdc, bvydc, bvzdc
  integer :: ix, iy, iz, ii, nnx,nny
  integer :: ix0, iy0, iz0
  real(long) :: rz0,r1,r2, xc,yc,zc,x1,x2,y1,y2,dx,dy
  real(long) :: px(1:12),py(1:12)
!
  x1 = -5.0d0*dis
  x2 = 5.0d0*dis
  y1 = -5.0d0*dis
  y2 = 5.0d0*dis
  dx = 0.5
  dy = 0.5
  nnx = (x2-x1)/dx
  nny = (y2-y1)/dy
!
  open(14,file='BBH_xy_vf.dat',status='unknown')
  zc = 0.0d0
  do iy = 1, nny+1
    yc = y1 + (iy-1)*dy
    write(14,*) ' '
    do ix = 1, nnx+1
      xc = x1 + (ix-1)*dx
      r1 = sqrt(xc**2 + yc**2 + zc**2)
      r2 = sqrt((xc-sepa)**2 + yc**2 + zc**2)
      if (r1.ge.rgin.and.r2.ge.rgin) then
        call interpo_gr2cgr_4th(bvxd,bvxdc,xc,yc,zc)
        call interpo_gr2cgr_4th(bvyd,bvydc,xc,yc,zc)
        call interpo_gr2cgr_4th(bvzd,bvzdc,xc,yc,zc)  
        if (dabs(bvxdc).le.small) bvxdc = 0.0d0
        if (dabs(bvydc).le.small) bvydc = 0.0d0
        if (dabs(bvzdc).le.small) bvzdc = 0.0d0
        write(14,'(20es14.6)') xc, yc, zc, bvxdc, bvydc, bvzdc
      else
        write(14,'(20es14.6)') xc, yc, zc, 0.0, 0.0, 0.0
      endif
    end do
  end do
  close(14)
!
  open(14,file='BBH_z1_vf.dat',status='unknown')
  zc = 0.5d0*rgin
  do iy = 1, nny+1
    yc = y1 + (iy-1)*dy
    write(14,*) ' '
    do ix = 1, nnx
      xc = x1 + (ix-1)*dx
      r1 = sqrt(xc**2 + yc**2 + zc**2)
      r2 = sqrt((xc-sepa)**2 + yc**2 + zc**2)
      if (r1.ge.rgin.and.r2.ge.rgin) then
        call interpo_gr2cgr_4th(bvxd,bvxdc,xc,yc,zc)
        call interpo_gr2cgr_4th(bvyd,bvydc,xc,yc,zc)
        call interpo_gr2cgr_4th(bvzd,bvzdc,xc,yc,zc)  
        if (dabs(bvxdc).le.small) bvxdc = 0.0d0
        if (dabs(bvydc).le.small) bvydc = 0.0d0
        if (dabs(bvzdc).le.small) bvzdc = 0.0d0
        write(14,'(20es14.6)') xc, yc, zc, bvxdc, bvydc, bvzdc
      else
        write(14,'(20es14.6)') xc, yc, zc, 0.0, 0.0, 0.0
      endif
    end do
  end do
  close(14)
!
  px(1) = -2.0d0*rgin
  py(1) = -2.0d0*rgin

  px(2) = 0.0d0
  py(2) = py(1)

  px(3) = 0.5d0*sepa
  py(3) = py(1)

  px(4) = sepa
  py(4) = py(1)

  px(5) = sepa + 2.0d0*rgin
  py(5) = py(1)

  px(6) = px(5)
  py(6) = 0.0d0

  px(7) = px(5)
  py(7) = 2.0d0*rgin

  px(8) = px(4)
  py(8) = py(7)

  px(9) = px(3)
  py(9) = py(7)

  px(10)= px(2)
  py(10)= py(7)

  px(11)= px(1)
  py(11)= py(7)

  px(12)= px(1)
  py(12)= 0.0d0
!
  open(14,file='BBH_plane_z0.dat',status='unknown')
  zc = 0.0d0
  do ii = 1, 12
    xc = px(ii)
    yc = py(ii)
    call interpo_gr2cgr_4th(bvxd,bvxdc,xc,yc,zc)
    call interpo_gr2cgr_4th(bvyd,bvydc,xc,yc,zc)
    call interpo_gr2cgr_4th(bvzd,bvzdc,xc,yc,zc)
    if (dabs(bvxdc).le.small) bvxdc = 0.0d0
    if (dabs(bvydc).le.small) bvydc = 0.0d0
    if (dabs(bvzdc).le.small) bvzdc = 0.0d0
    write(14,'(20es14.6)') xc, yc, zc, bvxdc, bvydc, bvzdc
    write(14,'(20es14.6)') xc, yc, zc, bvxdc, bvydc, 0.0d0
    write(14,'(20es14.6)') xc, yc, zc, 0.0d0, 0.0d0, bvzdc
  end do
 close(14)

 
!
end subroutine IO_output_cartesian_planes
