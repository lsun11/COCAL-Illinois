module grid_points_asymptotic_patch
  use phys_constant, only : long
  implicit none
!
  real(long), pointer :: ra(:,:,:), tha(:,:,:), phia(:,:,:)
  real(long), pointer :: hra(:,:,:), htha(:,:,:), hphia(:,:,:)
!
contains
!
subroutine allocate_grid_points_asymptotic_patch
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use make_array_3d
  implicit none
!
! -- ra, tha, phia, are extended coordinates ONLY in radial direction.
!
  call alloc_array3d(ra,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(tha,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(phia,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(hra,1,nrg,1,ntg,1,npg)
  call alloc_array3d(htha,1,nrg,1,ntg,1,npg)
  call alloc_array3d(hphia,1,nrg,1,ntg,1,npg)
!
end subroutine allocate_grid_points_asymptotic_patch
!
! --- Determin coordinates and grid number at the excision radius.
!
subroutine calc_grid_points_asymptotic_patch(impt_bin,impt_apt)
  use phys_constant, only : long, pi
  use coordinate_grav_r, only : rg, hrg
  use coordinate_grav_extended, only : rgex, hrgex
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only : nrg, ntg, npg, ntgxy, npgxzp, npgxzm, rgout
  use def_binary_parameter, only : dis
  use grid_parameter_binary_excision, only : ex_radius
  use trigonometry_grav_theta, only :  sinthg, costhg, hsinthg, hcosthg
  use trigonometry_grav_phi, only :  sinphig, cosphig, hsinphig, hcosphig
  implicit none
!
  real(long) :: small = 1.0d-10
  real(long) :: rrgg, sth, cth, sphi, cphi
  real(long) :: point_line_dis2, root_det, dis_bin
  real(long) :: x0, y0, z0, rg_exc_dis2, xxyy2, phi_rot
  integer    :: irg, itg, ipg, impt_bin, impt_apt
!
! --  coordinates of the aspymtptic grids seen from the centeral coordinates.
!
  call copy_def_binary_parameter_from_mpt(impt_bin)
  if (impt_bin.eq.1) dis_bin =   dis
  if (impt_bin.eq.2) dis_bin = - dis
  call copy_def_binary_parameter_from_mpt(impt_apt)
  phi_rot = 0.0d0
  if (impt_bin.eq.2) phi_rot = pi
!
! -- for impt=2 patch, rotate pi
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = -2, nrg + 2
        ra(irg,itg,ipg) = 0.0d0
        tha(irg,itg,ipg) = 0.0d0
        phia(irg,itg,ipg) = 0.0d0
!
        rrgg = rgex(irg)
        sth = sinthg(itg)
        cth = costhg(itg)
        sphi = sinphig(ipg)
        cphi = cosphig(ipg)
        x0 = rrgg*sth*cphi
        y0 = rrgg*sth*sphi
        z0 = rrgg*cth
!
        rg_exc_dis2 = (x0 + dis_bin)**2 + y0**2 + z0**2
        xxyy2       = (x0 + dis_bin)**2 + y0**2
!
        ra(irg,itg,ipg)  = sqrt(rg_exc_dis2)
        tha(irg,itg,ipg) = atan2(sqrt(xxyy2),z0)
        phia(irg,itg,ipg)= &
        &  dmod(2.0d0*pi+datan2(y0,x0+dis_bin)+phi_rot,2.0d0*pi)
!
      end do
    end do
  end do
!
! -- coordinates of the central half grid seen from the center of excision.
! -- for impt=2 patch, rotate pi
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        hra(irg,itg,ipg) = 0.0d0
        htha(irg,itg,ipg) = 0.0d0
        hphia(irg,itg,ipg) = 0.0d0
!
        rrgg = hrgex(irg)
        sth = hsinthg(itg)
        cth = hcosthg(itg)
        sphi = hsinphig(ipg)
        cphi = hcosphig(ipg)
        x0 = rrgg*sth*cphi
        y0 = rrgg*sth*sphi
        z0 = rrgg*cth
!
        rg_exc_dis2 = (x0 + dis_bin)**2 + y0**2 + z0**2
        xxyy2       = (x0 + dis_bin)**2 + y0**2
!
        hra(irg,itg,ipg)  = sqrt(rg_exc_dis2)
        htha(irg,itg,ipg) = atan2(sqrt(xxyy2),z0)
        hphia(irg,itg,ipg)= &
        &  dmod(2.0d0*pi+datan2(y0,x0+dis_bin)+phi_rot,2.0d0*pi)
!
      end do
    end do
  end do
!
end subroutine calc_grid_points_asymptotic_patch
end module grid_points_asymptotic_patch
