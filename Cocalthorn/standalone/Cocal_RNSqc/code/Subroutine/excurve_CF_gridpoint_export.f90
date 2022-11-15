subroutine excurve_CF_gridpoint_export(alph,bvxd,bvyd,bvzd,tfkxx,tfkxy,tfkxz,tfkyy,tfkyz,tfkzz)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
!  use def_metric, only : alph, bvxd, bvyd, bvzd
!  use def_metric_excurve_grid, only : tfkij_grid, tfkijkij_grid
  use interface_grgrad_4th_gridpoint
  implicit none
  integer :: ia, ib, info, irg, itg, ipg
  real(long) :: ainvh, diver, cdivbv, fa23, &
  &             dbvxdx, dbvxdy, dbvxdz, &
  &             dbvydx, dbvydy, dbvydz, &
  &             dbvzdx, dbvzdy, dbvzdz
  real(8), pointer :: alph(:,:,:), bvxd(:,:,:), bvyd(:,:,:), bvzd(:,:,:), tfkxx(:,:,:),  &
     &    tfkxy(:,:,:), tfkxz(:,:,:), tfkyy(:,:,:), tfkyz(:,:,:), tfkzz(:,:,:)
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
        call grgrad_4th_gridpoint(bvxd,dbvxdx,dbvxdy,dbvxdz,irg,itg,ipg)
        call grgrad_4th_gridpoint(bvyd,dbvydx,dbvydy,dbvydz,irg,itg,ipg)
        call grgrad_4th_gridpoint(bvzd,dbvzdx,dbvzdy,dbvzdz,irg,itg,ipg)
!
        ainvh = 0.5d0/alph(irg,itg,ipg)
        cdivbv = dbvxdx + dbvydy + dbvzdz
        diver = fa23*cdivbv
!
        tfkxx(irg,itg,ipg) = ainvh*(2.0d0*dbvxdx - diver) 
        tfkyy(irg,itg,ipg) = ainvh*(2.0d0*dbvydy - diver) 
        tfkzz(irg,itg,ipg) = ainvh*(2.0d0*dbvzdz - diver) 
        tfkxy(irg,itg,ipg) = ainvh*(dbvydx + dbvxdy)
        tfkxz(irg,itg,ipg) = ainvh*(dbvzdx + dbvxdz)
        tfkyz(irg,itg,ipg) = ainvh*(dbvzdy + dbvydz)
!        tfkij_grid(irg,itg,ipg,2,1) = tfkij_grid(irg,itg,ipg,1,2)
!        tfkij_grid(irg,itg,ipg,3,1) = tfkij_grid(irg,itg,ipg,1,3)
!        tfkij_grid(irg,itg,ipg,3,2) = tfkij_grid(irg,itg,ipg,2,3)
      end do
    end do
  end do
!
end subroutine excurve_CF_gridpoint_export
