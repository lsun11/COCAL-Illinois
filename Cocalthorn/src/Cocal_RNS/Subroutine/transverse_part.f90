subroutine transverse_part
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_transverse_part, only : Ftvx,      Ftvy,      Ftvz, &
  &                               Ftvx_grid, Ftvy_grid, Ftvz_grid
  use def_vector_x, only : vec_xg, hvec_xg
  use interface_grgrad_midpoint_type0
  use interface_grgrad_4th_gridpoint
  implicit none
  real(8) :: dhxxdx,dhxxdy,dhxxdz,dhxydx,dhxydy,dhxydz,dhxzdx,dhxzdy,dhxzdz, &
  &          dhyxdx,dhyxdy,dhyxdz,dhyydx,dhyydy,dhyydz,dhyzdx,dhyzdy,dhyzdz, &
  &          dhzxdx,dhzxdy,dhzxdz,dhzydx,dhzydy,dhzydz,dhzzdx,dhzzdy,dhzzdz
  real(8) :: x, y, z, dgamma(3)
  integer :: ipg, itg, irg
!
! --- Compute Transeverse part of spatial metric. 
!     whose value is assigned on the mid points. 
!     Should be used with routines to impose 
! --- gauge condition onto F^a = zD_b h^ab
!
! --- midpoints
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
! --- Derivatives of h_ab.
!
        call grgrad_midpoint_type0(hxxu,dhxxdx,dhxxdy,dhxxdz,irg,itg,ipg)
        call grgrad_midpoint_type0(hxyu,dhxydx,dhxydy,dhxydz,irg,itg,ipg)
        call grgrad_midpoint_type0(hxzu,dhxzdx,dhxzdy,dhxzdz,irg,itg,ipg)
        call grgrad_midpoint_type0(hyyu,dhyydx,dhyydy,dhyydz,irg,itg,ipg)
        call grgrad_midpoint_type0(hyzu,dhyzdx,dhyzdy,dhyzdz,irg,itg,ipg)
        call grgrad_midpoint_type0(hzzu,dhzzdx,dhzzdy,dhzzdz,irg,itg,ipg)
        dhyxdx = dhxydx ; dhyxdy = dhxydy ; dhyxdz = dhxydz
        dhzxdx = dhxzdx ; dhzxdy = dhxzdy ; dhzxdz = dhxzdz
        dhzydx = dhyzdx ; dhzydy = dhyzdy ; dhzydz = dhyzdz
!
        Ftvx(irg,itg,ipg) = dhxxdx + dhxydy + dhxzdz
        Ftvy(irg,itg,ipg) = dhyxdx + dhyydy + dhyzdz
        Ftvz(irg,itg,ipg) = dhzxdx + dhzydy + dhzzdz
!
      end do
    end do
  end do
!
! --- gridpoints
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
! --- Derivatives of h_ab.
!
        call grgrad_4th_gridpoint(hxxu,dhxxdx,dhxxdy,dhxxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hxyu,dhxydx,dhxydy,dhxydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hxzu,dhxzdx,dhxzdy,dhxzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hyyu,dhyydx,dhyydy,dhyydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hyzu,dhyzdx,dhyzdy,dhyzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(hzzu,dhzzdx,dhzzdy,dhzzdz,irg,itg,ipg)
        dhyxdx = dhxydx ; dhyxdy = dhxydy ; dhyxdz = dhxydz
        dhzxdx = dhxzdx ; dhzxdy = dhxzdy ; dhzxdz = dhxzdz
        dhzydx = dhyzdx ; dhzydy = dhyzdy ; dhzydz = dhyzdz
!
        Ftvx_grid(irg,itg,ipg) = dhxxdx + dhxydy + dhxzdz
        Ftvy_grid(irg,itg,ipg) = dhyxdx + dhyydy + dhyzdz
        Ftvz_grid(irg,itg,ipg) = dhzxdx + dhzydy + dhzzdz
!
      end do
    end do
  end do
!
!irg = 2; itg = ntg/2-2; ipg = 2
!x = hvec_xg(irg,itg,ipg,1) 
!y = hvec_xg(irg,itg,ipg,2)
!z = hvec_xg(irg,itg,ipg,3)
!write(6,'(1p,3e16.8)') x, y, z
!write(6,'(1p,3e16.8)') Ftvx(irg,itg,ipg), Ftvy(irg,itg,ipg), &
!&                      Ftvz(irg,itg,ipg)
!call kerr_schild_transverse_part(x,y,z,dgamma)
!write(6,'(1p,3e16.8)') dgamma(1), dgamma(2), dgamma(3)
!x = vec_xg(irg,itg,ipg,1) 
!y = vec_xg(irg,itg,ipg,2)
!z = vec_xg(irg,itg,ipg,3)
!write(6,'(1p,3e16.8)') x, y, z
!write(6,'(1p,3e16.8)') Ftvx_grid(irg,itg,ipg), Ftvy_grid(irg,itg,ipg), &
!&                      Ftvz_grid(irg,itg,ipg)
!call kerr_schild_transverse_part(x,y,z,dgamma)
!write(6,'(1p,3e16.8)') dgamma(1), dgamma(2), dgamma(3)
!stop
end subroutine transverse_part
