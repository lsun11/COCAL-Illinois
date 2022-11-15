subroutine excurve_CF_gridpoint_bhex
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : alph, bvxd, bvyd, bvzd
  use def_metric_excurve_grid, only : tfkij_grid, tfkijkij_grid, trk_grid
  use interface_grgrad_4th_gridpoint_bhex
  implicit none
  integer :: ia, ib, info, irg, itg, ipg
  real(long) :: ainvh, diver, cdivbv, fa23, &
  &             dbvxdx, dbvxdy, dbvxdz, &
  &             dbvydx, dbvydy, dbvydz, &
  &             dbvzdx, dbvzdy, dbvzdz
!
! --- Compute extringic curvature.  
! --- Whose value is assigned on the grid points. 
!
  fa23 = 2.0d0/3.0d0
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
!
        call grgrad_4th_gridpoint_bhex(bvxd,dbvxdx,dbvxdy,dbvxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(bvyd,dbvydx,dbvydy,dbvydz,irg,itg,ipg)
        call grgrad_4th_gridpoint_bhex(bvzd,dbvzdx,dbvzdy,dbvzdz,irg,itg,ipg)
!
        ainvh = 0.5d0/alph(irg,itg,ipg)
        cdivbv = dbvxdx + dbvydy + dbvzdz
        diver = fa23*cdivbv
!
        tfkij_grid(irg,itg,ipg,1,1) = ainvh*(2.0d0*dbvxdx - diver) 
        tfkij_grid(irg,itg,ipg,2,2) = ainvh*(2.0d0*dbvydy - diver) 
        tfkij_grid(irg,itg,ipg,3,3) = ainvh*(2.0d0*dbvzdz - diver) 
        tfkij_grid(irg,itg,ipg,1,2) = ainvh*(dbvydx + dbvxdy)
        tfkij_grid(irg,itg,ipg,1,3) = ainvh*(dbvzdx + dbvxdz)
        tfkij_grid(irg,itg,ipg,2,3) = ainvh*(dbvzdy + dbvydz)
        tfkij_grid(irg,itg,ipg,2,1) = tfkij_grid(irg,itg,ipg,1,2)
        tfkij_grid(irg,itg,ipg,3,1) = tfkij_grid(irg,itg,ipg,1,3)
        tfkij_grid(irg,itg,ipg,3,2) = tfkij_grid(irg,itg,ipg,2,3)
!
        tfkijkij_grid(irg,itg,ipg) = 0.0d0
        trk_grid(irg,itg,ipg) = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            tfkijkij_grid(irg,itg,ipg) = tfkijkij_grid(irg,itg,ipg) &
            &    +tfkij_grid(irg,itg,ipg,ia,ib)*tfkij_grid(irg,itg,ipg,ia,ib)
          end do
        end do
!
        if (tfkijkij_grid(irg,itg,ipg) /= 0.) info = 1
!
      end do
    end do
  end do
  if (info /= 1) write(6,*) ' ### Warning K_ij = 0 *** '
!
end subroutine excurve_CF_gridpoint_bhex
