subroutine calc_redblue_shift
  use phys_constant, only  : long, pi
  use coordinate_grav_r, only : rg
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfpolp, ntfeq, ntfxy, npfxzp, npfyzp
  use trigonometry_grav_theta, only : sinthg, costhg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd
  use def_matter, only : rs
  use def_matter_parameter, only : ome, ber, radi
  use def_quantities, only : zrb_xp_plus, zrb_xp_minus, &
  &                          zrb_yp_plus, zrb_yp_minus, &
  &                          zrb_zp_plus, zrb_zp_minus
  use interface_interpo_radial1p_grav
  implicit none
  integer    :: ir, it, ip, ii
  real(long) :: rv, xx, yy, zz, vphixu, vphiyu, vphizu
  real(long) :: psiw, alphw, bvxdw, bvydw, bvzdw
  real(long) :: gtt, gtp, gpp, hhgc, utgc, bplus, bminus
  real(long) :: kdott_plus, kdotu_plus, kdott_minus, kdotu_minus
!
! Semi-major axis along x
! Semi-minor axis along y
!
  do ii = 1, 2
    if (ii.eq.1) then; it = ntfeq; ip = npfxzp; end if
    if (ii.eq.2) then; it = ntfeq; ip = npfyzp; end if
!
    rv = rg(nrf)*rs(it,ip)
    xx = rv*sinthg(it)*cosphig(ip)
    yy = rv*sinthg(it)*sinphig(ip)
    zz = rv*costhg(it)
    vphixu = - yy
    vphiyu =   xx
    vphizu =   0.0d0
    hhgc  = 1.0d0
    utgc  = hhgc/ber
    call interpo_radial1p_grav(psi,psiw,rv,it,ip) 
    call interpo_radial1p_grav(alph,alphw,rv,it,ip) 
    call interpo_radial1p_grav(bvxd,bvxdw,rv,it,ip) 
    call interpo_radial1p_grav(bvyd,bvydw,rv,it,ip) 
    call interpo_radial1p_grav(bvzd,bvzdw,rv,it,ip) 
!
    gtt = - alphw**2 + psiw**4*(bvxdw**2 + bvydw**2 + bvzdw**2)
    gtp =   psiw**4*(bvxdw*vphixu + bvydw*vphiyu + bvzdw*vphizu)
    gpp =   psiw**4*(vphixu**2 + vphiyu**2 + vphizu**2)
    bplus = (- gtp + sqrt(gtp**2 - gtt*gpp))/gpp
    bminus = (- gtp - sqrt(gtp**2 - gtt*gpp))/gpp
!
    kdott_plus = gtt + bplus*gtp
    kdotu_plus = utgc*(gtt + (bplus + ome)*gtp + bplus*ome*gpp)
    kdott_minus = gtt + bminus*gtp
    kdotu_minus = utgc*(gtt + (bminus + ome)*gtp + bminus*ome*gpp)

    if (ii.eq.1) then
      zrb_xp_plus = kdotu_plus/kdott_plus - 1.0d0
      zrb_xp_minus = kdotu_minus/kdott_minus - 1.0d0
    end if
    if (ii.eq.2) then
      zrb_yp_plus = kdotu_plus/kdott_plus - 1.0d0
      zrb_yp_minus = kdotu_minus/kdott_minus - 1.0d0
    end if
!
  end do
!
  hhgc  = 1.0d0
  utgc  = hhgc/ber
  zrb_zp_plus = utgc - 1.0d0
  zrb_zp_minus = utgc - 1.0d0
!
end subroutine calc_redblue_shift
