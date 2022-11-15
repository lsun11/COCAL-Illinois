subroutine outer_boundary_d_psi(sou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  implicit none
  real(long), pointer :: sou_surf(:,:)
  integer    :: itg, ipg
! Reset boundary condition for BH
  do ipg = 1, npg
    do itg = 1, ntg
      sou_surf(itg,ipg) = 1.0d0
    end do
  end do
!
end subroutine outer_boundary_d_psi
