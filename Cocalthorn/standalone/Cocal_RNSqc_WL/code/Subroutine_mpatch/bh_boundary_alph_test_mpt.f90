subroutine bh_boundary_alph_test_mpt(impt,char_bc,sou_surf,dsou_surf)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use trigonometry_grav_theta, only : hsinthg
  use trigonometry_grav_phi,   only : hcosphig
  use def_binary_parameter,    only : sepa
  use grid_parameter_mpt,      only : grid_param_real_
  implicit none
  character(len=1), intent(in) :: char_bc
  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long) :: st, cp, rad1, rad2, dr2dr1, bhmass
  real(long) :: psi_tmp, dpsidr_tmp, alps_tmp, dalpsdr_tmp
  integer    :: itg, ipg
  integer    :: impt, impt_ex
  real(long) :: bhmass_ex, rgin_ex
!
  if (impt.eq.1) impt_ex = 2
  if (impt.eq.2) impt_ex = 1
  rgin_ex = grid_param_real_(1,impt_ex)
  bhmass_ex = 2.0d0*rgin_ex
  bhmass = 2.0d0*rgin
!
  do ipg = 1, npg
    do itg = 1, ntg
      st = hsinthg(itg)
      cp = hcosphig(ipg)
      rad1 = rgin
      rad2 = sqrt(rad1**2 - 2.0d0*rad1*sepa*st*cp + sepa**2)
      dr2dr1 = (rad1 - sepa*st*cp)/rad2
      psi_tmp  = 1.0d0 + 0.5d0*bhmass/rad1 + 0.5d0*bhmass_ex/rad2
      alps_tmp = 1.0d0 - 0.5d0*bhmass/rad1 - 0.5d0*bhmass_ex/rad2
      dpsidr_tmp = - 0.5d0*bhmass/rad1**2 &
      &            - 0.5d0*bhmass_ex/rad2**2*dr2dr1
      dalpsdr_tmp= + 0.5d0*bhmass/rad1**2 &
      &            + 0.5d0*bhmass_ex/rad2**2*dr2dr1
      if (char_bc.eq.'d') then
        sou_surf(itg,ipg) = alps_tmp/psi_tmp
      end if
      if (char_bc.eq.'n') then
        dsou_surf(itg,ipg) = dalpsdr_tmp/psi_tmp &
        &                  - alps_tmp*dpsidr_tmp/psi_tmp**2
      end if
    end do
  end do
!
end subroutine bh_boundary_alph_test_mpt
