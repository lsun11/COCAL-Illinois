subroutine interpolation_metric_bh(fnc,fncca)
  use phys_constant, only : long
  use coordinate_grav_r, only : rg
  use grid_parameter_cartesian, only : nx, ny, nz
  use coordinate_grav_xyz, only : x, y, z
  use interface_modules_cartesian, ignore_me => interpolation_metric_bh
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long), pointer :: fncca(:,:,:)
  real(long) :: xc, yc, zc, cfn, rc, r0
  integer :: ix, iy, iz
!
  r0 = rg(0)

  do iz = 1, nz
    zc = z(iz)
    do iy = 1, ny
      yc = y(iy)
      do ix = 1, nx
        xc = x(ix)
        rc = sqrt(xc**2 + yc**2 + zc**2)

        if ( rc>1.1d0*r0 ) then
          call interpo_gr2cgr_4th(fnc,cfn,xc,yc,zc)
          fncca(ix,iy,iz) = cfn
        else
          fncca(ix,iy,iz) = 0.0d0
        end if
      end do
    end do
  end do
end subroutine interpolation_metric_bh
