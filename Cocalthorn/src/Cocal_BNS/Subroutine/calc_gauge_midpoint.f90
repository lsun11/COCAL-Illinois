subroutine calc_gauge_midpoint
  use phys_constant, only : pi, long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : trk, ditrk
  use def_metric_excurve_grid, only : trk_grid
  use def_metric_dihiju
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use make_array_3d
  implicit none
!
  real(long), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(long) :: trkgc, Fxgc, Fygc, Fzgc
  integer :: ipg, irg, itg
!
  call alloc_array3d(dfdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdz,1,nrg,1,ntg,1,npg)
!
!
  call grgrad_midpoint(trk_grid, dfdx,dfdy,dfdz)
  ditrk(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  ditrk(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  ditrk(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(trkgc,trk_grid,   irg,itg,ipg)
        call interpo_linear_type0(Fxgc, dihixu_grid,irg,itg,ipg)
        call interpo_linear_type0(Fygc, dihiyu_grid,irg,itg,ipg)
        call interpo_linear_type0(Fzgc, dihizu_grid,irg,itg,ipg)

        trk(irg,itg,ipg)    = trkgc
        dihixu(irg,itg,ipg) = Fxgc
        dihiyu(irg,itg,ipg) = Fygc
        dihizu(irg,itg,ipg) = Fzgc
      end do
    end do
  end do
!
  deallocate(dfdx)
  deallocate(dfdy)
  deallocate(dfdz)
!
end subroutine calc_gauge_midpoint

