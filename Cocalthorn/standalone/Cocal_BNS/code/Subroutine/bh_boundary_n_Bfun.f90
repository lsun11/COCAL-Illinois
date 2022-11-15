subroutine bh_boundary_n_Bfun(dsou_surf,potx,poty,potz)
  use phys_constant, only  : long
  use grid_parameter, only : nrg, ntg, npg, rgin
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi,   only : hsinphig, hcosphig
  use def_binary_parameter,    only : sepa
  use def_metric
  implicit none
  real(long), pointer :: dsou_surf(:,:), potx(:,:,:), poty(:,:,:), potz(:,:,:)
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, tcm, pcm, xcm, ycm, zcm
  real(long) :: fx,fy,gx,gy
  integer    :: itg, ipg
!
  do ipg = 1, npg
    do itg = 1, ntg
      st = hsinthg(itg)
      ct = hcosthg(itg)
      sp = hsinphig(ipg)
      cp = hcosphig(ipg)
! Coordinates wrt the center of black hole
      xa = rgin*st*cp
      ya = rgin*st*sp
      za = rgin*ct
! Coordinates wrt the CM
      xcm = xa - 0.5d0*sepa
      ycm = ya
      zcm = za
 
      gx = 0.25d0*(potx(0,itg,ipg)   + potx(0,itg-1,ipg)  &
                 + potx(0,itg,ipg-1) + potx(0,itg-1,ipg-1) )
      gy = 0.25d0*(poty(0,itg,ipg)   + poty(0,itg-1,ipg)  &
                 + poty(0,itg,ipg-1) + poty(0,itg-1,ipg-1) )

      fx = 0.25d0*(bvxd(0,itg,ipg)   + bvxd(0,itg-1,ipg)  &
                 + bvxd(0,itg,ipg-1) + bvxd(0,itg-1,ipg-1) )
      fy = 0.25d0*(bvyd(0,itg,ipg)   + bvyd(0,itg-1,ipg)  &
                 + bvyd(0,itg,ipg-1) + bvyd(0,itg-1,ipg-1) )

      dsou_surf(itg,ipg) = 0.0d0
!      dsou_surf(itg,ipg) = 2.0d0*0.08*(-ycm*st*cp + xcm*st*sp)
!      dsou_surf(itg,ipg) = 4.0d0*(gx-fx)*st*cp + 4.0d0*(gy-fy)*st*sp
    end do
  end do
!
end subroutine bh_boundary_n_Bfun
