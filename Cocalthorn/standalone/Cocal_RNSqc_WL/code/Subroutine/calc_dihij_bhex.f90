subroutine calc_dihij_bhex
  use phys_constant, only :  long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_metric_dihiju 
  use def_formulation, only : chgra
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use interface_grgrad1g_midpoint
  use interface_grgrad_4th_gridpoint_bhex
  use make_array_3d
  implicit none
  real(long) :: grad1(3), dhu(3,3,3)
  real(long) :: dhyxdx, dhyxdy, dhyxdz, dhzxdx, dhzxdy, dhzxdz, &
     &          dhzydx, dhzydy, dhzydz, &
     &          dhxxdx, dhxxdy, dhxxdz, dhxydx, dhxydy, dhxydz, &
     &          dhxzdx, dhxzdy, dhxzdz, dhyydx, dhyydy, dhyydz, &
     &          dhyzdx, dhyzdy, dhyzdz, dhzzdx, dhzzdy, dhzzdz
  integer :: ipg, itg, irg, ia, ib, ic, id
!

  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call grgrad1g_midpoint(hxxu,grad1,irg,itg,ipg)
        dhu(1,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hxyu,grad1,irg,itg,ipg)
        dhu(1,2,1:3) = grad1(1:3)
        dhu(2,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hxzu,grad1,irg,itg,ipg)
        dhu(1,3,1:3) = grad1(1:3)
        dhu(3,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hyyu,grad1,irg,itg,ipg)
        dhu(2,2,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hyzu,grad1,irg,itg,ipg)
        dhu(2,3,1:3) = grad1(1:3)
        dhu(3,2,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hzzu,grad1,irg,itg,ipg)
        dhu(3,3,1:3) = grad1(1:3)
!
        dihixu(irg,itg,ipg) = dhu(1,1,1) + dhu(1,2,2) + dhu(1,3,3)
        dihiyu(irg,itg,ipg) = dhu(2,1,1) + dhu(2,2,2) + dhu(2,3,3)
        dihizu(irg,itg,ipg) = dhu(3,1,1) + dhu(3,2,2) + dhu(3,3,3)
      end do
    end do
  end do
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_4th_gridpoint_bhex(hxxu,dhxxdx,dhxxdy,dhxxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hxyu,dhxydx,dhxydy,dhxydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hxzu,dhxzdx,dhxzdy,dhxzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hyyu,dhyydx,dhyydy,dhyydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hyzu,dhyzdx,dhyzdy,dhyzdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(hzzu,dhzzdx,dhzzdy,dhzzdz,irg,itg,ipg)
!
        dihixu_grid(irg,itg,ipg) = dhxxdx + dhxydy + dhxzdz
        dihiyu_grid(irg,itg,ipg) = dhxydx + dhyydy + dhyzdz
        dihizu_grid(irg,itg,ipg) = dhxzdx + dhyzdy + dhzzdz
      end do
    end do
  end do
!
!
end subroutine calc_dihij_bhex
