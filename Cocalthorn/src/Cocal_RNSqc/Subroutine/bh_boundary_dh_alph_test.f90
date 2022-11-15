subroutine bh_boundary_dh_alph_test(sou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use trigonometry_grav_theta, only : hsinthg
  use trigonometry_grav_phi,   only : hcosphig
  use def_binary_parameter,    only : sepa
  implicit none
!  character(len=2), intent(in) :: char_bc
!  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long), pointer :: sou_surf(:,:)
  real(long) :: st, cp, rad1, rad2, dr2dr1, bhmass
  real(long) :: psisurf, alphsurf, dpsidr1surf
  integer    :: itg, ipg
!
! The factor 0.8 is to make alpha positive
  bhmass = 2.0d0*rgin*0.8d0


  do ipg = 1, npg
    do itg = 1, ntg
      st = hsinthg(itg)
      cp = hcosphig(ipg)
      rad1 = rgin
      rad2 = sqrt(rad1**2 - 2.0d0*rad1*sepa*st*cp + sepa**2)
!     dr2dr1 = (rad1 - sepa*st*cp)/rad2
! 
      psisurf = 1.0d0 + 0.5d0*bhmass/rad1 + 0.5d0*bhmass/rad2
!     dpsidr1surf = -0.5d0*bhmass/rad1**2 - 0.5d0*bhmass/rad2**2*dr2dr1
!     alphsurf = (1.0d0 - 0.5d0*bhmass/rad1 - 0.5d0*bhmass/rad2)/psisurf
      sou_surf(itg,ipg) = (1.0d0-0.5d0*bhmass/rad1-0.5d0*bhmass/rad2)/psisurf 
    end do
  end do
!
end subroutine bh_boundary_dh_alph_test
