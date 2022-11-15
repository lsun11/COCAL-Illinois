subroutine bh_boundary_1bh_nh_psi_test(dsou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use trigonometry_grav_theta, only : hsinthg
  use trigonometry_grav_phi,   only : hcosphig
  use def_binary_parameter,    only : sepa
  implicit none
!  character(len=2), intent(in) :: char_bc
!  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long), pointer :: dsou_surf(:,:)
   real(long) :: st, cp, rad1, bhmass
  integer    :: itg, ipg
!
  bhmass = 2.0d0*rgin
  do ipg = 1, npg
    do itg = 1, ntg
      st = hsinthg(itg)
      cp = hcosphig(ipg)
      rad1 = rgin
      dsou_surf(itg,ipg) = - 0.5d0*bhmass/rad1**2
    end do
  end do
!
end subroutine bh_boundary_1bh_nh_psi_test
