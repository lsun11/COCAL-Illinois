subroutine outer_boundary_d_psi_schw(sou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin, rgout
  use trigonometry_grav_theta, only : hsinthg
  use trigonometry_grav_phi,   only : hcosphig
  use def_binary_parameter,    only : dis
  use def_quantities,          only : admmass
  implicit none
  real(long), pointer :: sou_surf(:,:)
  real(long) :: st, cp, rad1, ri, bhmass
  integer    :: itg, ipg
! Reset boundary condition for BH
  bhmass = admmass
  do ipg = 1, npg
    do itg = 1, ntg
      st = hsinthg(itg)
      cp = hcosphig(ipg)
      rad1 = rgout
      ri   = sqrt(rad1**2 - 2.0d0*rad1*dis*st*cp + dis**2)
      sou_surf(itg,ipg) = 1.0d0 + 0.5d0*bhmass/ri
    end do
  end do
!
end subroutine outer_boundary_d_psi_schw
