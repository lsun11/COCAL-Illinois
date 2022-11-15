subroutine bh_boundary_d_potz(sou_surf,Bfun,dBfundz)
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, rgin
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi,   only : hsinphig, hcosphig
  use def_binary_parameter,    only : sepa
  implicit none
  real(long), pointer :: sou_surf(:,:), Bfun(:,:,:), dBfundz(:,:,:)
  real(long) :: dBdx,dBdy,dBdz, bf
  real(long) :: st, ct, sp, cp, xa,ya,za, rcm2, xycm2, rcm, tcm, pcm
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
      rcm2 = (xa-0.5d0*sepa)**2 + ya**2 + za**2
      xycm2= (xa-0.5d0*sepa)**2 + ya**2

      rcm = sqrt(rcm2)
      tcm = atan2(sqrt(xycm2),za)
      pcm = dmod(2.0d0*pi+datan2(ya,xa-0.5d0*sepa),2.0d0*pi)

      dBdz = 0.25d0*(dBfundz(0,itg,ipg)   + dBfundz(0,itg-1,ipg)  &
                   + dBfundz(0,itg,ipg-1) + dBfundz(0,itg-1,ipg-1) )

      bf = 0.25d0*(Bfun(0,itg,ipg)   + Bfun(0,itg-1,ipg)  &
                 + Bfun(0,itg,ipg-1) + Bfun(0,itg-1,ipg-1) )

!      sou_surf(itg,ipg) = 0.25d0*dBdz/bf
!      sou_surf(itg,ipg) = 0.25d0*dBdz
      sou_surf(itg,ipg) = 0.0d0
    end do
  end do
!
end subroutine bh_boundary_d_potz
