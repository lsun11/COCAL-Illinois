subroutine dadbscalar_type3(fnc,d2f,irg,itg,ipg)
  use phys_constant, only : long
  use interface_grgrad_4th_gridpoint
  use interface_grgrad_2nd_gridpoint
  implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long) :: d2f(1:3,1:3)
  real(long) :: dfdx, dfdy, dfdz
  real(long) :: dfncdx(0:1,0:1,0:1), dfncdy(0:1,0:1,0:1), dfncdz(0:1,0:1,0:1)
  real(long) :: d2fdxdx, d2fdxdy, d2fdxdz, &
              & d2fdydx, d2fdydy, d2fdydz, &
              & d2fdzdx, d2fdzdy, d2fdzdz
  integer :: irg, itg, ipg
  integer :: ip, ir, it, irg8, itg8, ipg8
!
! --- Compute the DaDb fnc of a function in 4th+2nd order.
! --- The DaDb fnc is evaluated at the mid points.
!
  do ip = 0, 1
    ipg8 = ipg - 1 + ip
    do it = 0, 1
      itg8 = itg - 1 + it
      do ir = 0, 1
        irg8 = irg - 1 + ir
!
!        if (irg.le.5) then 
!          call grgrad_2nd_gridpoint(fnc,dfdx,dfdy,dfdz,irg8,itg8,ipg8)
!        else
          call grgrad_4th_gridpoint(fnc,dfdx,dfdy,dfdz,irg8,itg8,ipg8)
!        end if
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
  d2f(1,1) = d2fdxdx
  d2f(1,2) =(d2fdxdy + d2fdydx)*0.5d0
  d2f(1,3) =(d2fdxdz + d2fdzdx)*0.5d0
  d2f(2,1) =(d2fdxdy + d2fdydx)*0.5d0
  d2f(2,2) = d2fdydy
  d2f(2,3) =(d2fdydz + d2fdzdy)*0.5d0
  d2f(3,1) =(d2fdxdz + d2fdzdx)*0.5d0
  d2f(3,2) =(d2fdydz + d2fdzdy)*0.5d0
  d2f(3,3) = d2fdzdz
!
end subroutine dadbscalar_type3
