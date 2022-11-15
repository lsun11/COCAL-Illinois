include '.././Module/phys_constant.f90'
include '.././Module/def_matter.f90'
include '.././Module/grid_parameter.f90'
include '.././Module/coordinate_grav_r.f90'
include '.././Module/coordinate_grav_phi.f90'
include '.././Module/coordinate_grav_theta.f90'
include '.././Module/coordinate_grav_extended.f90'
include '.././Module/make_array_1d.f90'
include '.././Module/make_array_2d.f90'
include '.././Module/make_array_3d.f90'
include '.././Module/make_array_4d.f90'
include '.././Module/make_array_5d.f90'
include '.././Module/trigonometry_grav_theta.f90'
include '.././Module/trigonometry_grav_phi.f90'
!
module interface_modules
  use phys_constant, only : long
  implicit none
  interface 
    subroutine grgrad_4th(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
      real(8), pointer :: fnc(:,:,:)
      real(8), intent(out) :: dfdx, dfdy, dfdz
      integer :: irg, itg, ipg
    end subroutine grgrad_4th
  end interface
end module interface_modules
!
include '.././Subroutine/grgrad_4th.f90'
include '.././Function/dfdx_4th.f90'
!##############################################
!
program main
implicit none
  call coordinate_patch_kit_grav
  call fncdiff
end program main
subroutine coordinate_patch_kit_grav
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use make_array_3d
  implicit none
! call subroutines. the order is important.
  call read_parameter
  call grid_r
  call grid_theta
  call trig_grav_theta
  call grid_phi
  call trig_grav_phi
  call grid_extended
end subroutine coordinate_patch_kit_grav
!
subroutine fncdiff
  use grid_parameter
  use coordinate_grav_r
  use coordinate_grav_phi
  use coordinate_grav_theta
  use coordinate_grav_extended
  use trigonometry_grav_theta
  use trigonometry_grav_phi
  use make_array_3d
  use interface_modules
  implicit none
  real(8), pointer :: fnc(:,:,:) 
  real(8) :: x, y, z, dfdx, dfdy, dfdz
  real(8) :: df1, df2, df3, e1, e2, e3, error, small = 1.0d-14
  integer :: irg, itg, ipg
!
  call alloc_array3d(fnc,0,nrg,0,ntg,0,npg)
  do irg = 0, 30
  do itg = 0, ntg
  do ipg = 0, npg
  x = rg(irg)*sinthg(itg)*cosphig(ipg)
  y = rg(irg)*sinthg(itg)*sinphig(ipg)
  z = rg(irg)*costhg(itg)
  fnc(irg,itg,ipg) = ((x**4 + 1)*y + y**3)*dcos(z)
  end do
  end do
  end do
!
  do ipg = 0, npg
  do itg = 0, ntg
  do irg = 0, 5
!
    x = rg(irg)*sinthg(itg)*cosphig(ipg)
    y = rg(irg)*sinthg(itg)*sinphig(ipg)
    z = rg(irg)*costhg(itg)
!
    call grgrad_4th(fnc,dfdx,dfdy,dfdz,irg,itg,ipg)
!
      write(6,*) irg,itg,ipg
      write(6,20) x,y,z
      write(6,20) dfdx,dfdy,dfdz
!
      df1 = 4*x**3*y*Cos(z)
      df2 = (1 + x**4 + 3*y**2)*Cos(z)
      df3 = -(((1 + x**4)*y + y**3)*Sin(z))
!
      e1 = 0.0d0
      e2 = 0.0d0
      e3 = 0.0d0
      error = 0.0d0
      if (dabs(df1)+dabs(dfdx).ge.small) & 
      e1 = dabs(df1-dfdx)/(dabs(df1)+dabs(dfdx))
      if (dabs(df2)+dabs(dfdy).ge.small) &
      e2 = dabs(df2-dfdy)/(dabs(df2)+dabs(dfdz))
      if (dabs(df3)+dabs(dfdz).ge.small) &
      e3 = dabs(df3-dfdz)/(dabs(df3)+dabs(dfdy))
      error = dmax1(e1,e2,e3)
      write(6,20) df1,df2,df3
      write(6,20) e1,e2,e3,error
  end do
  end do
  end do
!
  20 format(1p,6e12.4)
end subroutine fncdiff
