subroutine source_ang_mom_WL(soug,souf)
  use phys_constant, only :  long
  use grid_parameter, only :  nrg, ntg, npg, nrf, ntf, npf, ntgeq
  use coordinate_grav_r, only :   rg, hrg
  use trigonometry_grav_theta, only :  sinthg, costhg, hsinthg, hcosthg
  use trigonometry_grav_phi, only :  sinphig, cosphig, hsinphig, hcosphig
  use def_matter, only :  emd, rs, utf, omef
  use def_matter_parameter
  use def_metric, only : psi, alph, bvxd, bvyd, bvzd, tfkij, tfkijkij
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, &
  &                          hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use def_SEM_tensor, only : jmd
  use def_metric_on_SFC_CF
  use def_metric_on_SFC_WL
  use def_dvphi
  use make_array_3d
  use interface_interpo_gr2fl
  use interface_interpo_linear_type0
  use interface_grgrad1g_midpoint
  use def_vector_phi, only : hvec_phig, vec_phif
  implicit none
  real(long), pointer ::  soug(:,:,:),  souf(:,:,:)
  integer     ::   ir,it,ip
  real(long)  ::   emdw, alphw, psiw, rhow, prew, hhw, utw, omew, ene
  real(long)  ::   rjjx, rjjy, rjjz, rjjphi
  real(long)  ::   zfac, small = 1.0d-15
  real(long)  ::   vphif(1:3)
  real(long)  ::   otermx, otermy, otermz, bvxdfw, bvydfw, bvzdfw
  real(long)  ::   hhxxdf, hhxydf, hhxzdf, hhyxdf, hhyydf, hhyzdf, &
  &                hhzxdf, hhzydf, hhzzdf
  real(long)  ::   hxxuc, hxyuc, hxzuc, hyyuc, hyzuc, hzzuc
  real(long)  ::   gamu(3,3)
  real(long)  ::   aijvpp, aij(1:3,1:3), vphig(1:3), grad1(1:3), &
  &                dhu(1:3,1:3,1:3), psigc
  integer     ::   ipg, itg, irg, ia, ib
!
  dphiu(1:3,1:3) = 0.0d0
  dphiu(1,2) =-1.0d0
  dphiu(2,1) = 1.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(hxxuc,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hxyuc,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hxzuc,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hyyuc,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hyzuc,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hzzuc,hzzu,irg,itg,ipg)
        gamu(1,1) = hxxuc + 1.0d0
        gamu(1,2) = hxyuc
        gamu(1,3) = hxzuc
        gamu(2,2) = hyyuc + 1.0d0
        gamu(2,3) = hyzuc
        gamu(3,3) = hzzuc + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        aij(1:3,1:3) = tfkij(irg,itg,ipg,1:3,1:3)
        vphig(1) = hvec_phig(irg,itg,ipg,1)
        vphig(2) = hvec_phig(irg,itg,ipg,2)
        vphig(3) = hvec_phig(irg,itg,ipg,3)
        call grgrad1g_midpoint(hxxu,grad1,irg,itg,ipg)
        dhu(1,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hxyu,grad1,irg,itg,ipg)
        dhu(1,2,1:3) = grad1(1:3)
        dhu(2,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hxzu,grad1,irg,itg,ipg)
        dhu(1,3,1:3) = grad1(1:3)
        dhu(3,1,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hyyu,grad1,irg,itg,ipg)
        dhu(2,2,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hyzu,grad1,irg,itg,ipg)
        dhu(2,3,1:3) = grad1(1:3)
        dhu(3,2,1:3) = grad1(1:3)
        call grgrad1g_midpoint(hzzu,grad1,irg,itg,ipg)
        dhu(3,3,1:3) = grad1(1:3)
        aijvpp = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            aijvpp = aijvpp - 0.5d0*( aij(ia,ib)*vphig(1)*dhu(ia,ib,1)   &
                   &                + aij(ia,ib)*vphig(2)*dhu(ia,ib,2)   &
                   &                + aij(ia,ib)*vphig(3)*dhu(ia,ib,3) ) &
                   &                + aij(ia,ib)*gamu(ia,1)*dphiu(ib,1)  &
                   &                + aij(ia,ib)*gamu(ia,2)*dphiu(ib,2)  &
                   &                + aij(ia,ib)*gamu(ia,3)*dphiu(ib,3)
          end do
        end do
        soug(irg,itg,ipg) = aijvpp*psigc**6
      end do
    end do
  end do
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
!!        emdw = emd(ir,it,ip)
!!        if (emdw <= small) emdw = small
!!        psiw = psif(ir,it,ip)
!!        alphw = alphf(ir,it,ip)
!!        bvxdfw = bvxdf(ir,it,ip)
!!        bvydfw = bvydf(ir,it,ip)
!!        bvzdfw = bvzdf(ir,it,ip)
!!        call peos_q2hprho(emdw, hhw, prew, rhow, ene)
!!!
!!        utw  = utf(ir,it,ip)
!!        omew = omef(ir,it,ip)
!!!
!!        vphif(1) = vec_phif(ir,it,ip,1)
!!        vphif(2) = vec_phif(ir,it,ip,2)
!!        vphif(3) = vec_phif(ir,it,ip,3)
!!!
!!        hhxxdf = hxxdf(ir,it,ip)
!!        hhxydf = hxydf(ir,it,ip)
!!        hhxzdf = hxzdf(ir,it,ip)
!!        hhyydf = hyydf(ir,it,ip)
!!        hhyzdf = hyzdf(ir,it,ip)
!!        hhzzdf = hzzdf(ir,it,ip)
!!        hhyxdf = hhxydf
!!        hhzxdf = hhxzdf
!!        hhzydf = hhyzdf
!!!
!!        otermx = bvxdfw + omew*vphif(1) &
!!        &      +   hhxxdf*omew*vphif(1) &
!!        &      +   hhxydf*omew*vphif(2) &
!!        &      +   hhxzdf*omew*vphif(3)
!!        otermy = bvydfw + omew*vphif(2) &
!!        &      +   hhyxdf*omew*vphif(1) &
!!        &      +   hhyydf*omew*vphif(2) &
!!        &      +   hhyzdf*omew*vphif(3)
!!        otermz = bvzdfw + omew*vphif(3) &
!!        &      +   hhzxdf*omew*vphif(1) &
!!        &      +   hhzydf*omew*vphif(2) &
!!        &      +   hhzzdf*omew*vphif(3)
!!!
!!        zfac = 1.0d0
!!        if (emdw <= small) zfac = 0.0d0
!!        rjjx = hhw*rhow*alphw*utw**2*psiw**4*otermx
!!        rjjy = hhw*rhow*alphw*utw**2*psiw**4*otermy
!!        rjjz = hhw*rhow*alphw*utw**2*psiw**4*otermz
!!!
!!        rjjphi = rjjx*vphif(1) + rjjy*vphif(2) + rjjz*vphif(3)
!!        souf(ir,it,ip) = rjjphi*psiw**6*zfac
!
        psiw = psif(ir,it,ip)
        vphif(1) = vec_phif(ir,it,ip,1)
        vphif(2) = vec_phif(ir,it,ip,2)
        vphif(3) = vec_phif(ir,it,ip,3)
        rjjx = jmd(ir,it,ip,1)
        rjjy = jmd(ir,it,ip,2)
        rjjz = jmd(ir,it,ip,3)
!
        rjjphi = rjjx*vphif(1) + rjjy*vphif(2) + rjjz*vphif(3)
        souf(ir,it,ip) = rjjphi*psiw**6
!
      end do
    end do
  end do
!
end subroutine source_ang_mom_WL
