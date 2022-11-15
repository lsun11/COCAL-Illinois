subroutine grd2phi_midpoint_type0(fnc,deriv,irg,itg,ipg)
  use phys_constant, only : long
  use coordinate_grav_phi, only : hphig
  use coordinate_grav_extended, only : ipgex_phi, phigex
  implicit none
  real(long), external :: d2fdx2_2nd
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: deriv
  real(long) :: pv
  real(long) :: r4(4), th4(4), phi4(4), fp4(4)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ip0, ipg0, ii, jj, kk
!
! --- Compute the gradient of a function in 4th order.
! --- The gradient is evaluated on the grid points.
!
! --- r, theta, phi derivatives.
!
  ip0 = ipg-2
  deriv = 0.0d0
!
  do jj = 1, 2
    irgex = irg - 1 + jj
    do kk = 1, 2
      itgex = itg - 1 + kk
      do ii = 1, 4
        ipg0 = ip0 + ii - 1
        phi4(ii) = phigex(ipg0)
!
        ipgex = ipgex_phi(ipg0)
        fp4(ii) = fnc(irgex,itgex,ipgex)
      end do
!
! --- To cartesian component.
!
      pv = hphig(ipg)
      deriv = deriv + d2fdx2_2nd(phi4,fp4,pv)
    end do
  end do
!
  deriv = 0.25d0*deriv
!
end subroutine grd2phi_midpoint_type0
