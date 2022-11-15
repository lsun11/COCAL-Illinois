subroutine IO_output_cartesian_contour_potential_BBH_CF_v1
  use phys_constant, only : long
  use def_metric, only : psi,alph
  use grid_parameter_cartesian, only  : nx, ny, nz, nx_mid
  use grid_parameter_binary_excision, only : ex_rgmid, ex_radius 
  use coordinate_grav_xyz, only  : x, y, z
  use def_binary_parameter,    only : sepa, dis
  use grid_parameter, only : rgin
  use interface_modules_cartesian
  implicit none
  real(long) :: small = 1.0d-20
  real(long) :: psic, alphc, bvxdc, bvydc, bvzdc
  integer :: ix, iy, iz, ii, nnx,nny
  integer :: ix0, iy0, iz0, npoints, counter
  real(long) :: rz0,r1,r2, xc,yc,zc,x1,x2,y1,y2,dx,dy
!
  x1 = -5.0d0*dis
  x2 = 5.0d0*dis
  y1 = -5.0d0*dis
  y2 = 5.0d0*dis
  nnx= 1500
  nny= 1500
  dx = (x2-x1)/nnx
  dy = (y2-y1)/nny

  npoints = 0
  yc = 0.0d0
! Change zc to get contour plots on another z=const plane
  zc = 0.0d0
  if(dabs(zc).lt.rgin) then
    rz0 = sqrt(rgin**2 - zc**2)
  else 
    rz0 =0.0d0
  endif
  do ix = 1, nnx+1
    xc = x1 + (ix-1)*dx
    r1 = dabs(xc)
    r2 = dabs(xc-sepa)
    if (r1.ge.rz0.and.r2.ge.rz0) then
      npoints = npoints + 1
    endif
  end do

  open(14,file='BBH_contour_xy_v1.dat',status='unknown')
  zc = 0.0d0
  do iy = 1, nny+1
    yc = y1 + (iy-1)*dy
    counter = 0
    write(14,*) ' '
    do ix = 1, nnx+1
      xc = x1 + (ix-1)*dx
      r1 = sqrt(xc**2 + yc**2 + zc**2)
      r2 = sqrt((xc-sepa)**2 + yc**2 + zc**2)
      if (r1.ge.rgin.and.r2.ge.rgin) then
        if (r2 <= ex_radius*1.2d0) then
          call interpo_gr2cgr_4th(psi,psic,-xc+ex_rgmid,-yc,zc)
          call interpo_gr2cgr_4th(alph,alphc,-xc+ex_rgmid,-yc,zc)
        else
          call interpo_gr2cgr_4th(psi,psic,xc,yc,zc)
          call interpo_gr2cgr_4th(alph,alphc,xc,yc,zc)
        endif
          counter = counter + 1 
          write(14,'(20es14.6)') xc, yc, psic, alphc
          if (counter==npoints) exit 
      endif
    end do
  end do
  close(14)
!
end subroutine IO_output_cartesian_contour_potential_BBH_CF_v1
