module grid_points_binary_in_asympto
  use phys_constant, only : long
  implicit none
!
  real(long), pointer :: rb_a(:,:,:), thb_a(:,:,:), phib_a(:,:,:)
  real(long), pointer :: hrb_a(:,:,:), hthb_a(:,:,:), hphib_a(:,:,:)
!
contains
!
subroutine allocate_grid_points_binary_in_asympto
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use make_array_3d
  implicit none
!
! -- rb_a, thb_a, phib_a, are extended coordinates ONLY in radial direction.
!
  call alloc_array3d(rb_a,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(thb_a,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(phib_a,-2,nrg+2,0,ntg,0,npg)
  call alloc_array3d(hrb_a,1,nrg,1,ntg,1,npg)
  call alloc_array3d(hthb_a,1,nrg,1,ntg,1,npg)
  call alloc_array3d(hphib_a,1,nrg,1,ntg,1,npg)
!
end subroutine allocate_grid_points_binary_in_asympto
!
! --- Determin coordinates and grid number at the excision radius.
!
subroutine calc_grid_points_binary_in_asympto(impt_bin,impt_apt)
  use phys_constant, only : long, pi
  use coordinate_grav_r, only : rg, hrg
  use coordinate_grav_extended, only : rgex, hrgex
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use grid_parameter, only : nrg, ntg, npg, ntgxy, npgxzp, npgxzm, rgin
  use def_binary_parameter, only : dis
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
! -- coordinates of the binary grids seen from the asymptotic coordinates.
! -- for impt=2 patch, rotate pi
!
  dis_bin = dis
  phi_rot = 0.0d0
  if (impt_bin.eq.2) phi_rot = pi
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = -2, nrg + 2
        rb_a(irg,itg,ipg) = 0.0d0
        thb_a(irg,itg,ipg) = 0.0d0
        phib_a(irg,itg,ipg) = 0.0d0
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
        rg_exc_dis2 = (x0 - dis_bin)**2 + y0**2 + z0**2
        xxyy2       = (x0 - dis_bin)**2 + y0**2
!
        if (rg_exc_dis2.ge.rgin) then 
          rb_a(irg,itg,ipg)  = sqrt(rg_exc_dis2)
          thb_a(irg,itg,ipg) = atan2(sqrt(xxyy2),z0)
          phib_a(irg,itg,ipg)= &
          &  dmod(2.0d0*pi+datan2(y0,x0+dis_bin)+phi_rot,2.0d0*pi)
        end if
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
        hrb_a(irg,itg,ipg) = 0.0d0
        hthb_a(irg,itg,ipg) = 0.0d0
        hphib_a(irg,itg,ipg) = 0.0d0
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
        rg_exc_dis2 = (x0 - dis_bin)**2 + y0**2 + z0**2
        xxyy2       = (x0 - dis_bin)**2 + y0**2
!
        if (rg_exc_dis2.lt.rgin) then 
          hrb_a(irg,itg,ipg)  = sqrt(rg_exc_dis2)
          hthb_a(irg,itg,ipg) = atan2(sqrt(xxyy2),z0)
          hphib_a(irg,itg,ipg)= &
          &  dmod(2.0d0*pi+datan2(y0,x0+dis_bin)+phi_rot,2.0d0*pi)
        end if
!
      end do
    end do
  end do
!
end subroutine calc_grid_points_binary_in_asympto
end module grid_points_binary_in_asympto
