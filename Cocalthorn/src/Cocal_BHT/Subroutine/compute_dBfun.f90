subroutine compute_dBfun(Bfun,dBfundx,dBfundy,dBfundz)
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg
  use interface_grgrad_4th_gridpoint
  implicit none
  real(long) :: dBdx, dBdy, dBdz 
  real(long), pointer :: Bfun(:,:,:), dBfundx(:,:,:), dBfundy(:,:,:), dBfundz(:,:,:)
  integer :: irg, itg, ipg
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        call grgrad_4th_gridpoint(Bfun,dBdx,dBdy,dBdz,irg,itg,ipg)
        dBfundx(irg,itg,ipg) = dBdx
        dBfundy(irg,itg,ipg) = dBdy
        dBfundz(irg,itg,ipg) = dBdz
      end do
    end do
  end do
!
end subroutine compute_dBfun
