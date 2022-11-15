subroutine outer_boundary_d_potx(sou_surf,Bfun,dBfundx)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  implicit none
  real(long), pointer :: sou_surf(:,:), Bfun(:,:,:), dBfundx(:,:,:)
  real(long) :: dBdx,dBdy,dBdz, bf
  integer    :: itg, ipg
! Reset boundary condition for BH
  do ipg = 1, npg
    do itg = 1, ntg
      dBdx = 0.25d0*(dBfundx(nrg,itg,ipg)   + dBfundx(nrg,itg-1,ipg)  &
                   + dBfundx(nrg,itg,ipg-1) + dBfundx(nrg,itg-1,ipg-1) )

      bf = 0.25d0*(Bfun(0,itg,ipg)   + Bfun(0,itg-1,ipg)  &
                 + Bfun(0,itg,ipg-1) + Bfun(0,itg-1,ipg-1) )
!
!      sou_surf(itg,ipg) = 0.0d0 + 0.25d0*dBdx/bf
!      sou_surf(itg,ipg) = 0.0d0 + 0.25d0*dBdx
      sou_surf(itg,ipg) = 0.0d0 
    end do
  end do
!
end subroutine outer_boundary_d_potx
