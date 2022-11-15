subroutine calc_shift2rotshift
  use grid_parameter, only : nrg, ntg, npg
  use def_matter, only : omeg
  use def_matter_parameter, only : ber, radi
  use def_metric, only : bvxu, bvyu, bvzu, bvxd, bvyd, bvzd
  use coordinate_grav_r, only : rg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use trigonometry_grav_theta, only : sinthg
  use def_metric_rotshift, only : ovxu, ovyu, ovzu, ovxd, ovyd, ovzd
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_vector_phi, only : vec_phig
  implicit none
  real(8) :: vphixg, vphiyg, vphizg, bvxgc, bvygc, bvzgc, &
  &          gmxxd, gmxyd, gmxzd, gmyxd, gmyyd, gmyzd, &
  &          gmzxd, gmzyd, gmzzd, omegc
  integer :: ipg, itg, irg
!
! --- shift to rotating shift.  
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        vphixg = vec_phig(irg,itg,ipg,1)
        vphiyg = vec_phig(irg,itg,ipg,2)
        vphizg = vec_phig(irg,itg,ipg,3)
        bvxgc = bvxu(irg,itg,ipg)
        bvygc = bvyu(irg,itg,ipg)
        bvzgc = bvzu(irg,itg,ipg)
        omegc = omeg(irg,itg,ipg)
        ovxu(irg,itg,ipg) = bvxgc + omegc*vphixg
        ovyu(irg,itg,ipg) = bvygc + omegc*vphiyg
        ovzu(irg,itg,ipg) = bvzgc + omegc*vphizg
      end do
    end do
  end do
!
! --- Lowering index of rotating shift.  
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        gmxxd = 1.0d0 + hxxd(irg,itg,ipg)
        gmxyd =         hxyd(irg,itg,ipg)
        gmxzd =         hxzd(irg,itg,ipg)
        gmyxd =         hxyd(irg,itg,ipg)
        gmyyd = 1.0d0 + hyyd(irg,itg,ipg)
        gmyzd =         hyzd(irg,itg,ipg)
        gmzxd =         hxzd(irg,itg,ipg)
        gmzyd =         hyzd(irg,itg,ipg)
        gmzzd = 1.0d0 + hzzd(irg,itg,ipg)
!
        vphixg = vec_phig(irg,itg,ipg,1)
        vphiyg = vec_phig(irg,itg,ipg,2)
        vphizg = vec_phig(irg,itg,ipg,3)
        bvxgc = bvxd(irg,itg,ipg)
        bvygc = bvyd(irg,itg,ipg)
        bvzgc = bvzd(irg,itg,ipg)
        omegc = omeg(irg,itg,ipg)
!
        ovxd(irg,itg,ipg) = bvxgc + gmxxd*omegc*vphixg &
           &                      + gmxyd*omegc*vphiyg &
           &                      + gmxzd*omegc*vphizg
        ovyd(irg,itg,ipg) = bvygc + gmyxd*omegc*vphixg &
           &                      + gmyyd*omegc*vphiyg &
           &                      + gmyzd*omegc*vphizg
        ovzd(irg,itg,ipg) = bvzgc + gmzxd*omegc*vphixg &
           &                      + gmzyd*omegc*vphiyg &
           &                      + gmzzd*omegc*vphizg
      end do
    end do
  end do
!
end subroutine calc_shift2rotshift
!bv2ov
