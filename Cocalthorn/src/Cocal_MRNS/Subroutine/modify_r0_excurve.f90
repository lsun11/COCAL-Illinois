subroutine modify_r0_excurve
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only :  tfkij
  use def_metric_excurve_grid, only : tfkij_grid
  use coordinate_grav_r, only : rg
  implicit none
  real(long), external :: lagint_2nd
  real(long) :: x(2),y(2), v
  integer    :: irg, itg, ipg, ia, ib
!
  do ia = 1,3
    do ib = 1,3
      do ipg = 0,npg
        do itg = 0,ntg
          x(1) = rg(1)
          x(2) = rg(2)
          y(1) = tfkij_grid(1,itg,ipg,ia,ib)
          y(2) = tfkij_grid(2,itg,ipg,ia,ib)
          v = rg(0)      
          tfkij_grid(0,itg,ipg,ia,ib) = lagint_2nd(x,y,v)
        end do
      end do
    end do
  end do
!
end subroutine modify_r0_excurve

