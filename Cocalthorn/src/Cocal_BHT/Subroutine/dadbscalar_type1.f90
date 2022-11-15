subroutine dadbscalar_type1(fnc,d2fdxdx,d2fdxdy,d2fdxdz,d2fdydx,d2fdydy,d2fdydz,d2fdzdx,d2fdzdy,d2fdzdz,irg,itg,ipg)
  use phys_constant, only : long
  use interface_grgrad_2nd
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long) :: dfdx, dfdy, dfdz
  real(long) :: dfncdx(0:1,0:1,0:1), dfncdy(0:1,0:1,0:1), dfncdz(0:1,0:1,0:1)
  real(long) :: d2fdxdx, d2fdxdy, d2fdxdz, &
              & d2fdydx, d2fdydy, d2fdydz, &
              & d2fdzdx, d2fdzdy, d2fdzdz
  integer :: irg, itg, ipg
  integer :: irg0 , itg0 , ipg0, ii
  integer :: ip, ir, it
!  real(long), external :: dfdx_2nd
!
! --- Compute the gradient of a function in 3rd order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
! --- To cartesian component.
!
  do ip = 0, 1
    ipg0 = ipg - 1 + ip
    do it = 0, 1
      itg0 = itg - 1 + it
      do ir = 0, 1
        irg0 = irg - 1 + ir
        call grgrad_2nd(fnc,dfdx,dfdy,dfdz,irg0,itg0,ipg0)
!
        dfncdx(ir,it,ip) = dfdx
        dfncdy(ir,it,ip) = dfdy
        dfncdz(ir,it,ip) = dfdz
      end do
    end do
  end do
!
  call grgrad_type0(dfncdx,d2fdxdx,d2fdxdy,d2fdxdz,irg,itg,ipg)
  call grgrad_type0(dfncdy,d2fdydx,d2fdydy,d2fdydz,irg,itg,ipg)
  call grgrad_type0(dfncdz,d2fdzdx,d2fdzdy,d2fdzdz,irg,itg,ipg)
  d2fdxdy = (d2fdxdy + d2fdydx)*0.5d0
  d2fdydx =  d2fdxdy
  d2fdxdz = (d2fdxdz + d2fdzdx)*0.5d0
  d2fdzdx =  d2fdxdz
  d2fdydz = (d2fdydz + d2fdzdy)*0.5d0
  d2fdzdy =  d2fdydz
!
end subroutine dadbscalar_type1
