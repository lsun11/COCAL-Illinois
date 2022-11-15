subroutine compute_shift_v2(Bfun,potx,poty,potz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : bvxd, bvyd, bvzd
  use interface_grgrad_4th_gridpoint
  implicit none
  real(long) :: dBfundx, dBfundy, dBfundz 
  real(long), pointer :: Bfun(:,:,:), potx(:,:,:), poty(:,:,:), potz(:,:,:)
  integer :: irg, itg, ipg
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_4th_gridpoint(Bfun,dBfundx,dBfundy,dBfundz,irg,itg,ipg)
        bvxd(irg,itg,ipg) = potx(irg,itg,ipg) - 0.25d0*dBfundx
        bvyd(irg,itg,ipg) = poty(irg,itg,ipg) - 0.25d0*dBfundy
        bvzd(irg,itg,ipg) = potz(irg,itg,ipg) - 0.25d0*dBfundz
      end do
    end do
  end do
!
end subroutine compute_shift_v2
