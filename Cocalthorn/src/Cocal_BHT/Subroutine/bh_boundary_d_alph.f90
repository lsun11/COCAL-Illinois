subroutine bh_boundary_d_alph(sou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use def_binary_parameter,    only : sepa
  implicit none
  real(long), pointer :: sou_surf(:,:)
  integer    :: itg, ipg
!
  do ipg = 1, npg
    do itg = 1, ntg
      sou_surf(itg,ipg) = 1.0 
    end do
  end do
!
end subroutine bh_boundary_d_alph
