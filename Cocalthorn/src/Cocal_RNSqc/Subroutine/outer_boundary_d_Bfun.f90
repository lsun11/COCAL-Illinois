subroutine outer_boundary_d_Bfun(sou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  implicit none
  real(long), pointer :: sou_surf(:,:)
  integer    :: itg, ipg
!
  do ipg = 1, npg
    do itg = 1, ntg
!      sou_surf(itg,ipg) = 1.0d0
      sou_surf(itg,ipg) = 0.0d0
    end do
  end do
!
end subroutine outer_boundary_d_Bfun
