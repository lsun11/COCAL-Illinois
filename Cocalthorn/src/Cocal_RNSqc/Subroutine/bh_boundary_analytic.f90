subroutine bh_boundary_analytic(char_mp,irg,sou_surf,dsou_surf)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r,       only : rg
  use trigonometry_grav_theta, only : hsinthg,  hcosthg
  use trigonometry_grav_phi,   only : hsinphig, hcosphig
  implicit none
  real(long), pointer :: sou_surf(:,:), dsou_surf(:,:)
  real(long) :: psi, alph, bvu(3), bvd(3), hijd(3,3), hiju(3,3)
  real(long) :: x, y, z, sw(11)
  integer :: irg, itg, ipg
  character(len=4), intent(in) :: char_mp
!
! -- Kerr-Schild values at the black hole boundary and outer boundary.
! -- grid point for r, and mid point for theta and phi coordinates.
!
  dsou_surf(1:ntg,1:npg) = 0.0d0
  sw(1:11) = 0.0d0
  if (char_mp.eq.'psi ') sw(1) = 1.0d0
  if (char_mp.eq.'alph') sw(2) = 1.0d0
  if (char_mp.eq.'bvxd') sw(3) = 1.0d0
  if (char_mp.eq.'bvyd') sw(4) = 1.0d0
  if (char_mp.eq.'bvzd') sw(5) = 1.0d0
  if (char_mp.eq.'hxxd') sw(6) = 1.0d0
  if (char_mp.eq.'hxyd') sw(7) = 1.0d0
  if (char_mp.eq.'hxzd') sw(8) = 1.0d0
  if (char_mp.eq.'hyyd') sw(9) = 1.0d0
  if (char_mp.eq.'hyzd') sw(10)= 1.0d0
  if (char_mp.eq.'hzzd') sw(11)= 1.0d0
  do ipg = 1, npg
    do itg = 1, ntg
      x = rg(irg)*hsinthg(itg)*hcosphig(ipg)
      y = rg(irg)*hsinthg(itg)*hsinphig(ipg)
      z = rg(irg)*hcosthg(itg)
      call kerr_schild_metric_3plus1(x,y,z,psi,alph,bvu,bvd,hijd,hiju)
      sou_surf(itg,ipg) = sw(1)*psi + sw(2)*alph &
      &                 + sw(3)*bvd(1) + sw(4)*bvd(2) + sw(5)*bvd(3) &
      &                 + sw( 6)*hijd(1,1) + sw( 7)*hijd(1,2) &
      &                 + sw( 8)*hijd(1,3) + sw( 9)*hijd(2,2) &
      &                 + sw(10)*hijd(2,3) + sw(11)*hijd(3,3)
!!  if (char_mp.eq.'psi ') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'alph') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'bvxd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'bvyd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'bvzd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'hxxd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'hxyd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'hxzd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'hyyd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'hyzd') sou_surf(itg,ipg) = 1.0d0
!!  if (char_mp.eq.'hzzd') sou_surf(itg,ipg) = 1.0d0
    end do
  end do
!
end subroutine bh_boundary_analytic
