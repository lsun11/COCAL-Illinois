subroutine sourceterm_trfreeG_WL_bhex(souten)
  use grid_parameter, only : long, nrg, ntg, npg, nrf
  use phys_constant, only : pi
  use coordinate_grav_r, only : rg, hrg
  use def_matter_parameter, only : ome, ber, radi
  use def_metric, only : psi, alph, alps, tfkij, trk, &
  &                      bvxu, bvyu, bvzu, bvxd, bvyd, bvzd, alps2
  use def_metric_excurve_grid
  use def_CTT_decomposition, only : tfkij_grid_CTT
  use def_metric_rotshift, only : ovxu, ovyu, ovzu, &
  &                               ovxd, ovyd, ovzd
  use def_cristoffel, only : cri
  use def_ricci_tensor, only : rabnl, rab2d
  use def_shift_derivatives, only : cdbvxd, cdbvyd, cdbvzd, cdivbv
  use def_Lie_derivatives, only : rlpxx, rlpxy, rlpxz, rlpyy, rlpyz, rlpzz, &
  &                               rlbxx, rlbxy, rlbxz, rlbyy, rlbyz, rlbzz
  use def_Lie_derivatives_grid, only : rlbxx_grid, rlbxy_grid, rlbxz_grid, &
  &                                    rlbyy_grid, rlbyz_grid, rlbzz_grid, &
  &                                    rlpxx_grid, rlpxy_grid, rlpxz_grid, &
  &                                    rlpyy_grid, rlpyz_grid, rlpzz_grid
  use def_kerr_schild
  use def_bh_parameter, only : bh_bctype
  use def_formulation, only : swlp, swls, iswl, chgra, chope
  use def_cutsw, only : cutfac
  use def_dvphi, only : dphiu
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu, &
  &                          hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use trigonometry_grav_theta, only : hsinthg, hcosthg
  use trigonometry_grav_phi, only : hsinphig, hcosphig
  use make_array_3d
  use make_array_4d
  use make_array_5d
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use interface_grgrad1g_midpoint
  use interface_grdphi_midpoint_type0
  use interface_dadbscalar_type0
  use interface_dadbscalar_type3_bhex
  use interface_grd2phi_midpoint_type0
  implicit none
!
  real(8), pointer :: souten(:,:,:,:)
  real(8), pointer :: dfdx(:,:,:), dfdy(:,:,:), dfdz(:,:,:)
  real(8), pointer :: fnc0(:,:,:)
  real(8), pointer :: grada(:,:,:,:), gradp(:,:,:,:)
  real(8), pointer :: gradap(:,:,:,:), gradap2(:,:,:,:)
  real(8), pointer :: gradk(:,:,:,:,:), gradb(:,:,:,:,:)
  real(8) :: dabalps2(1:3,1:3), daalps2(3), &
  &          dapsi(3), daalph(3), daalps(3), cdbv(3,3)
  real(8) :: dbv(1:3,1:3), dov(1:3,1:3)=0.0d0
  real(8) :: daij(1:6,1:3)
  real(8) :: gamu(1:3,1:3), gamd(1:3,1:3), aij(1:3,1:3)
  real(8) :: rlieaij(3,3), r2dab(3,3), rnlab(3,3), rpsiab(3,3), liebaij(3,3), &
  &          aacabc(3,3), sab(3,3), ovd(3), ovu(3), rpddal(3,3), aijdb(3,3)
  real(8) :: heli(3,3), dphabd(3,3), dphabu(3,3), d2phab(3,3), &
  &          rlbgab(3,3), rlpgab(3,3), dliepg(6,3), dliebg(6,3), &
  &          rlplbg(3,3), rlblpg(3,3), rlblbg(3,3), palidp(3,3)
  real(8) :: ica(6), icb(6)
  real(8) :: tm1a(3,3),tm2a(3,3),tm3a(3,3),tm4a(3,3),tm5a(3,3)
!
  real(8) :: alpgc, pregc, psigc, &
  &          hhxxu, hhxyu, hhxzu, hhyxu, hhyyu, hhyzu, &
  &          hhzxu, hhzyu, hhzzu, &
  &          bvxdc, bvxgc, bvydc, bvygc, bvzdc, bvzgc, &
  &          san, san2, cutoff, &
  &          alpgcinv, ap2gc, ap2gcinv, c1ab, c2ab, c3ab, cdalps2, cdivbc, &
  &          d2phxx, d2phxy, d2phxz, d2phyx, d2phyy, d2phyz, &
  &          d2phzx, d2phzy, d2phzz, &
  &          dapda, dhdh, &
  &          dhxxddx, dhxxddy, dhxxddz, dhxxudx, dhxxudy, dhxxudz, &
  &          dhxyddx, dhxyddy, dhxyddz, dhxyudx, dhxyudy, dhxyudz, &
  &          dhxzddx, dhxzddy, dhxzddz, dhxzudx, dhxzudy, dhxzudz, &
  &          dhyxddx, dhyxddy, dhyxddz, dhyxudx, dhyxudy, dhyxudz, &
  &          dhyyddx, dhyyddy, dhyyddz, dhyyudx, dhyyudy, dhyyudz, &
  &          dhyzddx, dhyzddy, dhyzddz, dhyzudx, dhyzudy, dhyzudz, &
  &          dhzxddx, dhzxddy, dhzxddz, dhzxudx, dhzxudy, dhzxudz, &
  &          dhzyddx, dhzyddy, dhzyddz, dhzyudx, dhzyudy, dhzyudz, &
  &          dhzzddx, dhzzddy, dhzzddz, dhzzudx, dhzzudy, dhzzudz, &
  &          dpdbgxx, dpdbgxy, dpdbgxz, dpdbgyx, dpdbgyy, dpdbgyz, &
  &          dpdbgzx, dpdbgzy, dpdbgzz, &
  &          dphdph, dphxxd, dphxyd, dphxzd, dphxxu, dphxyu, dphxzu, &
  &          dphyxd, dphyyd, dphyzd, dphyxu, dphyyu, dphyzu, &
  &          dphzxd, dphzyd, dphzzd, dphzxu, dphzyu, dphzzu, &
  &          gadadp, gamdhdh, hhxxd, hhxyd, hhxzd, &
  &          ogdphdph, ovpadpsi, ovxgc, ovygc, ovzgc, bepadpsi, &
  &          ps43oal, ps4oal, ps4oal2, psigc4, psiinv, psiinv2, &
  &          rlielna, rlielnp, roku, sgam, sw1, sw2, sw3, &
  &          term1, term2, term3, term4, term5, tfaa, &
  &          tfaij, tfd2ph, tfheli, tflieaij, tfliep4aij, tfrnlab, tfrpddal, &
  &          trsab, traa, traij, trd2ph, trheli, trlieaij, trrnlab, trrpddal, &
  &          tfliebaij, trliebaij, trr2dab, tfr2dab,   &
  &          hhyxd, hhyyd, hhyzd, hhzxd, hhzyd, hhzzd, sou_daij,  &
  &          ene, grad(1:3), bda, bdp
  real(8) :: zerofac=0.0d0, xx,yy,zz,rr,rxy
  real(8) :: br,bth,bph, psi4, rr2, rxy2
  integer :: ipg, irg, itg, iph, ii, ia, ib, ic, irkij
!
  call alloc_array3d(dfdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfdz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(fnc0,0,nrg,0,ntg,0,npg)
  call alloc_array4d(grada,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array4d(gradp,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array4d(gradap,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array4d(gradap2,1,nrg,1,ntg,1,npg,1,3)
  call alloc_array5d(gradk,1,nrg,1,ntg,1,npg,1,6,1,3)
  call alloc_array5d(gradb,1,nrg,1,ntg,1,npg,1,3,1,3)
!
! --- source for hij is evaluated on grid points.
!
  do ii = 1, 3
! ### bvx, bvy, bvz are non-rotating shift!!
    if (ii == 1) call grgrad_midpoint(bvxu,dfdx,dfdy,dfdz)
    if (ii == 2) call grgrad_midpoint(bvyu,dfdx,dfdy,dfdz)
    if (ii == 3) call grgrad_midpoint(bvzu,dfdx,dfdy,dfdz)
    gradb(1:nrg,1:ntg,1:npg,ii,1) = dfdx(1:nrg,1:ntg,1:npg)
    gradb(1:nrg,1:ntg,1:npg,ii,2) = dfdy(1:nrg,1:ntg,1:npg)
    gradb(1:nrg,1:ntg,1:npg,ii,3) = dfdz(1:nrg,1:ntg,1:npg)
  end do
!
!
! Until irg=irkij take derivatives of tfkij to be that of Kerr Schild
  irkij = nrf  !10 close to converge but finally blows up,   !nrg  

  do ic = 1, 6
    ia = 1 + ic/4 + ic/6
    ib = ic - (ic/4)*2 - ic/6

!    if (bh_bctype.eq.'KS') then
!    if (chgra == 'k') then
      fnc0(0:nrg,0:ntg,0:npg) = tfkij_grid_CTT(0:nrg,0:ntg,0:npg,ia,ib)  
!    else
!      fnc0(0:nrg,0:ntg,0:npg) = tfkij_grid(0:nrg,0:ntg,0:npg,ia,ib)  
!    end if

!    fnc0(0:irkij,0:ntg,0:npg)     = tfkij_grid_ks(0:irkij,0:ntg,0:npg,ia,ib)
!    fnc0(irkij+1:nrg,0:ntg,0:npg) = tfkij_grid(irkij+1:nrg,0:ntg,0:npg,ia,ib)

    call grgrad_midpoint(fnc0,dfdx,dfdy,dfdz)
    gradk(1:nrg,1:ntg,1:npg,ic,1) = dfdx(1:nrg,1:ntg,1:npg)
    gradk(1:nrg,1:ntg,1:npg,ic,2) = dfdy(1:nrg,1:ntg,1:npg)
    gradk(1:nrg,1:ntg,1:npg,ic,3) = dfdz(1:nrg,1:ntg,1:npg)
  end do

!
  alps2(0:nrg,0:ntg,0:npg) = &
  &        alph(0:nrg,0:ntg,0:npg)*psi(0:nrg,0:ntg,0:npg)**2
!
  call grgrad_midpoint(psi,dfdx,dfdy,dfdz)
  gradp(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  gradp(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  gradp(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
  call grgrad_midpoint(alph,dfdx,dfdy,dfdz)
  grada(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  grada(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  grada(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
  call grgrad_midpoint(alps,dfdx,dfdy,dfdz)
  gradap(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  gradap(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  gradap(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
  call grgrad_midpoint(alps2,dfdx,dfdy,dfdz)
  gradap2(1:nrg,1:ntg,1:npg,1) = dfdx(1:nrg,1:ntg,1:npg)
  gradap2(1:nrg,1:ntg,1:npg,2) = dfdy(1:nrg,1:ntg,1:npg)
  gradap2(1:nrg,1:ntg,1:npg,3) = dfdz(1:nrg,1:ntg,1:npg)
!
  san = 1.0d0/3.0d0
  san2= 2.0d0/3.0d0
  roku= 1.0d0/6.0d0
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        cutoff = 1.0d0
        if (chgra == 'C'.or.chgra == 'H'.or.chgra == 'W') then
          if (rg(irg) > cutfac*pi/ome) cutoff = 0.0d0
        end if
!
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
        call interpo_linear_type0(ap2gc,alps2,irg,itg,ipg)
        psigc4 = psigc**4
        psiinv = 1.0d0/psigc
        psiinv2 = psiinv**2
        alpgcinv = 1.0d0/alpgc
        ap2gcinv = 1.0d0/ap2gc
        ps43oal = 4.0d0*psigc**3*alpgcinv
        ps4oal = psigc4*alpgcinv
        ps4oal2= psigc4*alpgcinv**2
!
        call interpo_linear_type0(hhxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hhxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hhxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hhyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hhyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hhzzu,hzzu,irg,itg,ipg)
        hhyxu = hhxyu
        hhzxu = hhxzu
        hhzyu = hhyzu
        gamu(1,1) = hhxxu + 1.0d0
        gamu(1,2) = hhxyu
        gamu(1,3) = hhxzu
        gamu(2,2) = hhyyu + 1.0d0
        gamu(2,3) = hhyzu
        gamu(3,3) = hhzzu + 1.0d0
        gamu(2,1) = gamu(1,2)
        gamu(3,1) = gamu(1,3)
        gamu(3,2) = gamu(2,3)
!
        call interpo_linear_type0(hhxxd,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hhxyd,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hhxzd,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hhyyd,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hhyzd,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hhzzd,hzzd,irg,itg,ipg)
        hhyxd = hhxyd
        hhzxd = hhxzd
        hhzyd = hhyzd
        gamd(1,1) = hhxxd + 1.0d0
        gamd(1,2) = hhxyd
        gamd(1,3) = hhxzd
        gamd(2,2) = hhyyd + 1.0d0
        gamd(2,3) = hhyzd
        gamd(3,3) = hhzzd + 1.0d0
        gamd(2,1) = gamd(1,2)
        gamd(3,1) = gamd(1,3)
        gamd(3,2) = gamd(2,3)
!
        call interpo_linear_type0(ovxgc,ovxd,irg,itg,ipg)
        call interpo_linear_type0(ovygc,ovyd,irg,itg,ipg)
        call interpo_linear_type0(ovzgc,ovzd,irg,itg,ipg)
        ovd(1) = ovxgc
        ovd(2) = ovygc
        ovd(3) = ovzgc
        call interpo_linear_type0(ovxgc,ovxu,irg,itg,ipg)
        call interpo_linear_type0(ovygc,ovyu,irg,itg,ipg)
        call interpo_linear_type0(ovzgc,ovzu,irg,itg,ipg)
        ovu(1) = ovxgc
        ovu(2) = ovygc
        ovu(3) = ovzgc
!
! ### bvx, bvy, bvz are non-rotating shift!!
        call interpo_linear_type0(bvxgc,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvygc,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzgc,bvzu,irg,itg,ipg)
!
        dapsi(1:3) = gradp(irg,itg,ipg,1:3)
        daalph(1:3) = grada(irg,itg,ipg,1:3)
        daalps(1:3) = gradap(irg,itg,ipg,1:3)
        daalps2(1:3) = gradap2(irg,itg,ipg,1:3)
        cdbv(1,1:3) = cdbvxd(irg,itg,ipg,1:3)
        cdbv(2,1:3) = cdbvyd(irg,itg,ipg,1:3)
        cdbv(3,1:3) = cdbvzd(irg,itg,ipg,1:3)
        do ia = 1, 3
          aij(ia,1:3) = tfkij(irg,itg,ipg,ia,1:3)
          dbv(ia,1:3) = gradb(irg,itg,ipg,ia,1:3)
        end do
        daij(1:6,1:3) = gradk(irg,itg,ipg,1:6,1:3)
!          if(itg.eq.ntg/2.and.ipg.eq.1) &
!          & write(6,*) 'daij', daij(1,1)
!
        cdivbc = cdivbv(irg,itg,ipg)
!
! --- DaDb(alpha psi^2)
!
!        call dadbscalar_type0(alps2,dabalps2,irg,itg,ipg)
        call dadbscalar_type3_bhex(alps2,dabalps2,irg,itg,ipg)
!
        gadadp=gamu(1,1)*daalph(1)*dapsi(1)+gamu(1,2)*daalph(1)*dapsi(2) &
        &     +gamu(1,3)*daalph(1)*dapsi(3)+gamu(2,1)*daalph(2)*dapsi(1) &
        &     +gamu(2,2)*daalph(2)*dapsi(2)+gamu(2,3)*daalph(2)*dapsi(3) &
        &     +gamu(3,1)*daalph(3)*dapsi(1)+gamu(3,2)*daalph(3)*dapsi(2) &
        &     +gamu(3,3)*daalph(3)*dapsi(3)
!
! --- Lie(psi^4 A_ab)
!
!        sw1 = swlp(iswl)
!        sw2 = 1.0d0 - sw1
!        sw3 = swls(iswl)
!
        bepadpsi=ps43oal*(bvxgc*dapsi(1)+bvygc*dapsi(2)+bvzgc*dapsi(3))

! ### rotating shift
!
        ovpadpsi=ps43oal*(ovxgc*dapsi(1)+ovygc*dapsi(2)+ovzgc*dapsi(3))
        dov(1:3,1:3) = dbv(1:3,1:3)
        dov(1,2) = dbv(1,2) + ome*dphiu(1,2)
        dov(2,1) = dbv(2,1) + ome*dphiu(2,1)
!
! --- For helical terms
!
        if (chgra == 'h'.or.chgra == 'c'.or.chgra == 'C' &
           &.or.chgra == 'H'.or.chgra == 'W') then
!
          rlbgab(1,1) = rlbxx(irg,itg,ipg)
          rlbgab(1,2) = rlbxy(irg,itg,ipg)
          rlbgab(1,3) = rlbxz(irg,itg,ipg)
          rlbgab(2,1) = rlbxy(irg,itg,ipg)
          rlbgab(2,2) = rlbyy(irg,itg,ipg)
          rlbgab(2,3) = rlbyz(irg,itg,ipg)
          rlbgab(3,1) = rlbxz(irg,itg,ipg)
          rlbgab(3,2) = rlbyz(irg,itg,ipg)
          rlbgab(3,3) = rlbzz(irg,itg,ipg)
          rlpgab(1,1) = rlpxx(irg,itg,ipg)
          rlpgab(1,2) = rlpxy(irg,itg,ipg)
          rlpgab(1,3) = rlpxz(irg,itg,ipg)
          rlpgab(2,1) = rlpxy(irg,itg,ipg)
          rlpgab(2,2) = rlpyy(irg,itg,ipg)
          rlpgab(2,3) = rlpyz(irg,itg,ipg)
          rlpgab(3,1) = rlpxz(irg,itg,ipg)
          rlpgab(3,2) = rlpyz(irg,itg,ipg)
          rlpgab(3,3) = rlpzz(irg,itg,ipg)
!
          call grdphi_midpoint_type0(hxxd,dphxxd,irg,itg,ipg)
          call grdphi_midpoint_type0(hxyd,dphxyd,irg,itg,ipg)
          call grdphi_midpoint_type0(hxzd,dphxzd,irg,itg,ipg)
          call grdphi_midpoint_type0(hyyd,dphyyd,irg,itg,ipg)
          call grdphi_midpoint_type0(hyzd,dphyzd,irg,itg,ipg)
          call grdphi_midpoint_type0(hzzd,dphzzd,irg,itg,ipg)
          dphyxd = dphxyd
          dphzxd = dphxzd
          dphzyd = dphyzd
          dphabd(1,1) = dphxxd
          dphabd(1,2) = dphxyd
          dphabd(1,3) = dphxzd
          dphabd(2,1) = dphyxd
          dphabd(2,2) = dphyyd
          dphabd(2,3) = dphyzd
          dphabd(3,1) = dphzxd
          dphabd(3,2) = dphzyd
          dphabd(3,3) = dphzzd
          call grdphi_midpoint_type0(hxxu,dphxxu,irg,itg,ipg)
          call grdphi_midpoint_type0(hxyu,dphxyu,irg,itg,ipg)
          call grdphi_midpoint_type0(hxzu,dphxzu,irg,itg,ipg)
          call grdphi_midpoint_type0(hyyu,dphyyu,irg,itg,ipg)
          call grdphi_midpoint_type0(hyzu,dphyzu,irg,itg,ipg)
          call grdphi_midpoint_type0(hzzu,dphzzu,irg,itg,ipg)
          dphyxu = dphxyu
          dphzxu = dphxzu
          dphzyu = dphyzu
          dphabu(1,1) = dphxxu
          dphabu(1,2) = dphxyu
          dphabu(1,3) = dphxzu
          dphabu(2,1) = dphyxu
          dphabu(2,2) = dphyyu
          dphabu(2,3) = dphyzu
          dphabu(3,1) = dphzxu
          dphabu(3,2) = dphzyu
          dphabu(3,3) = dphzzu
          call grd2phi_midpoint_type0(hxxd,d2phxx,irg,itg,ipg)
          call grd2phi_midpoint_type0(hxyd,d2phxy,irg,itg,ipg)
          call grd2phi_midpoint_type0(hxzd,d2phxz,irg,itg,ipg)
          call grd2phi_midpoint_type0(hyyd,d2phyy,irg,itg,ipg)
          call grd2phi_midpoint_type0(hyzd,d2phyz,irg,itg,ipg)
          call grd2phi_midpoint_type0(hzzd,d2phzz,irg,itg,ipg)
          d2phyx = d2phxy
          d2phzx = d2phxz
          d2phzy = d2phyz
          d2phab(1,1) = d2phxx
          d2phab(1,2) = d2phxy
          d2phab(1,3) = d2phxz
          d2phab(2,1) = d2phyx
          d2phab(2,2) = d2phyy
          d2phab(2,3) = d2phyz
          d2phab(3,1) = d2phzx
          d2phab(3,2) = d2phzy
          d2phab(3,3) = d2phzz
!
! --  \Lie_\phi T_{ab} 
! --   = \pa_\phi T_{ab} + T_{ac} \zD_b\phi^c + T_{cb} \zD_a\phi^c
!
          call grdphi_midpoint_type0(rlbxx_grid,dpdbgxx,irg,itg,ipg)
          call grdphi_midpoint_type0(rlbxy_grid,dpdbgxy,irg,itg,ipg)
          call grdphi_midpoint_type0(rlbxz_grid,dpdbgxz,irg,itg,ipg)
          call grdphi_midpoint_type0(rlbyy_grid,dpdbgyy,irg,itg,ipg)
          call grdphi_midpoint_type0(rlbyz_grid,dpdbgyz,irg,itg,ipg)
          call grdphi_midpoint_type0(rlbzz_grid,dpdbgzz,irg,itg,ipg)
          dpdbgyx = dpdbgxy
          dpdbgzx = dpdbgxz
          dpdbgzy = dpdbgyz
          rlplbg(1,1)=dpdbgxx+rlbgab(1,2)*dphiu(2,1)+rlbgab(2,1)*dphiu(2,1)
          rlplbg(1,2)=dpdbgxy+rlbgab(1,1)*dphiu(1,2)+rlbgab(2,2)*dphiu(2,1)
          rlplbg(1,3)=dpdbgxz+rlbgab(2,3)*dphiu(2,1)
          rlplbg(2,1)=dpdbgyx+rlbgab(2,2)*dphiu(2,1)+rlbgab(1,1)*dphiu(1,2)
          rlplbg(2,2)=dpdbgyy+rlbgab(2,1)*dphiu(1,2)+rlbgab(1,2)*dphiu(1,2)
          rlplbg(2,3)=dpdbgyz+rlbgab(1,3)*dphiu(1,2)
          rlplbg(3,1) = dpdbgzx
          rlplbg(3,2) = dpdbgzy
          rlplbg(3,3) = dpdbgzz
!
          rlielnp = 8.0d0*psiinv &
          &       *(ovu(1)* dapsi(1) + ovu(2)* dapsi(2) + ovu(3)* dapsi(3))
          rlielna = alpgcinv &
          &       *(ovu(1)*daalph(1) + ovu(2)*daalph(2) + ovu(3)*daalph(3))
!
          call grgrad1g_midpoint(rlbxx_grid,grad,irg,itg,ipg)
          dliebg(1,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlbxy_grid,grad,irg,itg,ipg)
          dliebg(2,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlbxz_grid,grad,irg,itg,ipg)
          dliebg(3,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlbyy_grid,grad,irg,itg,ipg)
          dliebg(4,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlbyz_grid,grad,irg,itg,ipg)
          dliebg(5,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlbzz_grid,grad,irg,itg,ipg)
          dliebg(6,1:3) = grad(1:3)
!
          call grgrad1g_midpoint(rlpxx_grid,grad,irg,itg,ipg)
          dliepg(1,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlpxy_grid,grad,irg,itg,ipg)
          dliepg(2,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlpxz_grid,grad,irg,itg,ipg)
          dliepg(3,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlpyy_grid,grad,irg,itg,ipg)
          dliepg(4,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlpyz_grid,grad,irg,itg,ipg)
          dliepg(5,1:3) = grad(1:3)
          call grgrad1g_midpoint(rlpzz_grid,grad,irg,itg,ipg)
          dliepg(6,1:3) = grad(1:3)
!
        end if
!
        do ic = 1, 6
!
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
!
! --- Lie_w A_{ab} 
!       = w^c \pa_c A_{ab} + A_{ac}\pa_b w^c + A_{cb}\pa_a w^c. 
!
          rlieaij(ia,ib) = &
          &     ovxgc*daij(ic,1) + ovygc*daij(ic,2) + ovzgc*daij(ic,3) &
          &   + aij(ia,1)*dov(1,ib)+aij(ia,2)*dov(2,ib)+aij(ia,3)*dov(3,ib) &
          &   + aij(1,ib)*dov(1,ia)+aij(2,ib)*dov(2,ia)+aij(3,ib)*dov(3,ia)

          liebaij(ia,ib) = &
          &     bvxgc*daij(ic,1) + bvygc*daij(ic,2) + bvzgc*daij(ic,3) &
          &   + aij(ia,1)*dbv(1,ib)+aij(ia,2)*dbv(2,ib)+aij(ia,3)*dbv(3,ib) &
          &   + aij(1,ib)*dbv(1,ia)+aij(2,ib)*dbv(2,ia)+aij(3,ib)*dbv(3,ia)

!         Following term is LieAij without the dAij terms
          aijdb(ia,ib) = &
          &   + aij(ia,1)*dbv(1,ib)+aij(ia,2)*dbv(2,ib)+aij(ia,3)*dbv(3,ib) &
          &   + aij(1,ib)*dbv(1,ia)+aij(2,ib)*dbv(2,ia)+aij(3,ib)*dbv(3,ia)



!          if(itg.eq.ntg/2.and.ipg.eq.1.and.ic.eq.1) then
!           write(6,*) 'rlieaij', rlieaij(1,1), aij(1,ib), aij(2,ib), aij(3,ib), dov(3,ia)
!          end if
!
! --- R_ab(nonlinear) 
!
          rnlab(ia,ib) = rabnl(irg,itg,ipg,ic)
!
! --- R_ab(Kerr_Schild) 
!
          r2dab(ia,ib) = rab2d(irg,itg,ipg,ic) 
!
!
! --- R^psi_ab - 1/alpha DaDb alpha
!
          c1ab = cri(irg,itg,ipg,1,ic)
          c2ab = cri(irg,itg,ipg,2,ic)
          c3ab = cri(irg,itg,ipg,3,ic)
!
          cdalps2 = c1ab*daalps2(1) + c2ab*daalps2(2) + c3ab*daalps2(3)
          dapda = 4.0d0*(daalps(ia)*dapsi(ib)+daalps(ib)*dapsi(ia))
          rpddal(ia,ib) = ap2gcinv*(- dabalps2(ia,ib) + cdalps2 + dapda)
!
! --- AacAbc
!
          aacabc(ia,ib) = gamu(1,1)*aij(ia,1)*aij(ib,1) &
          &             + gamu(1,2)*aij(ia,1)*aij(ib,2) &
          &             + gamu(1,3)*aij(ia,1)*aij(ib,3) &
          &             + gamu(2,1)*aij(ia,2)*aij(ib,1) &
          &             + gamu(2,2)*aij(ia,2)*aij(ib,2) &
          &             + gamu(2,3)*aij(ia,2)*aij(ib,3) &
          &             + gamu(3,1)*aij(ia,3)*aij(ib,1) &
          &             + gamu(3,2)*aij(ia,3)*aij(ib,2) &
          &             + gamu(3,3)*aij(ia,3)*aij(ib,3)
!
! --- Helical terms
!
          heli(ia,ib) = 0.0d0
          if (chgra == 'h'.or.chgra == 'H'.or.chgra == 'W') then
!
            rlblpg(ia,ib) = &
            &     bvxgc*dliepg(ic,1)+bvygc*dliepg(ic,2)+bvzgc*dliepg(ic,3) &
            &   + rlpgab(ia,1)*dov(1,ib) + rlpgab(ia,2)*dov(2,ib) &
            &   + rlpgab(ia,3)*dov(3,ib) + rlpgab(1,ib)*dov(1,ia) &
            &   + rlpgab(2,ib)*dov(2,ia) + rlpgab(3,ib)*dov(3,ia)
            rlblbg(ia,ib) = &
            &     bvxgc*dliebg(ic,1)+bvygc*dliebg(ic,2)+bvzgc*dliebg(ic,3) &
            &   + rlbgab(ia,1)*dov(1,ib) + rlbgab(ia,2)*dov(2,ib) &
            &   + rlbgab(ia,3)*dov(3,ib) + rlbgab(1,ib)*dov(1,ia) &
            &   + rlbgab(2,ib)*dov(2,ia) + rlbgab(3,ib)*dov(3,ia)
!
            palidp(ia,ib) = (dphabd(ia,1) + rlpgab(ia,1))*dphiu(1,ib) &
            &             + (dphabd(ia,2) + rlpgab(ia,2))*dphiu(2,ib) &
            &             + (dphabd(ia,3) + rlpgab(ia,3))*dphiu(3,ib) &
            &             + (dphabd(1,ib) + rlpgab(1,ib))*dphiu(1,ia) &
            &             + (dphabd(2,ib) + rlpgab(2,ib))*dphiu(2,ia) &
            &             + (dphabd(3,ib) + rlpgab(3,ib))*dphiu(3,ia)
!
            term1 = 0.5d0*(ps4oal2 - 1.0d0)*ome**2*d2phab(ia,ib)
            term2 = 0.5d0*ps4oal2*ome**2*palidp(ia,ib)
            term3 = 0.5d0*ps4oal2*ome*(rlplbg(ia,ib) + rlblpg(ia,ib))
            term4 = 0.5d0*ps4oal2*rlblbg(ia,ib)
            term5 = ps4oal*aij(ia,ib)*(rlielnp - rlielna)
!
! --  for print out
            if (irg == 2.and.itg == 1.and.ipg == 0) then
              write(6,'(1p,6e14.6)')term1,term2,term3,term4,term5
            end if
!
            tm1a(ia,ib) = term1
            tm2a(ia,ib) = term2
            tm3a(ia,ib) = term3
            tm4a(ia,ib) = term4
            tm5a(ia,ib) = term5
! --  
!
            if (chgra == 'H') then
              term2 = term2*cutoff
              term3 = term3*cutoff
            end if
!
            heli(ia,ib) = term1 + term2 + term3 + term4 + term5
!
          end if
!
        end do
!
! --- derivative for h_ab and h^ab.
!
        call grgrad1g_midpoint(hxxd,grad,irg,itg,ipg)
        dhxxddx = grad(1)
        dhxxddy = grad(2)
        dhxxddz = grad(3)
        call grgrad1g_midpoint(hxyd,grad,irg,itg,ipg)
        dhxyddx = grad(1)
        dhxyddy = grad(2)
        dhxyddz = grad(3)
        call grgrad1g_midpoint(hxzd,grad,irg,itg,ipg)
        dhxzddx = grad(1)
        dhxzddy = grad(2)
        dhxzddz = grad(3)
        call grgrad1g_midpoint(hyyd,grad,irg,itg,ipg)
        dhyyddx = grad(1)
        dhyyddy = grad(2)
        dhyyddz = grad(3)
        call grgrad1g_midpoint(hyzd,grad,irg,itg,ipg)
        dhyzddx = grad(1)
        dhyzddy = grad(2)
        dhyzddz = grad(3)
        call grgrad1g_midpoint(hzzd,grad,irg,itg,ipg)
        dhzzddx = grad(1)
        dhzzddy = grad(2)
        dhzzddz = grad(3)
!
        dhyxddx = dhxyddx
        dhyxddy = dhxyddy
        dhyxddz = dhxyddz
        dhzxddx = dhxzddx
        dhzxddy = dhxzddy
        dhzxddz = dhxzddz
        dhzyddx = dhyzddx
        dhzyddy = dhyzddy
        dhzyddz = dhyzddz
!
        call grgrad1g_midpoint(hxxu,grad,irg,itg,ipg)
        dhxxudx = grad(1)
        dhxxudy = grad(2)
        dhxxudz = grad(3)
        call grgrad1g_midpoint(hxyu,grad,irg,itg,ipg)
        dhxyudx = grad(1)
        dhxyudy = grad(2)
        dhxyudz = grad(3)
        call grgrad1g_midpoint(hxzu,grad,irg,itg,ipg)
        dhxzudx = grad(1)
        dhxzudy = grad(2)
        dhxzudz = grad(3)
        call grgrad1g_midpoint(hyyu,grad,irg,itg,ipg)
        dhyyudx = grad(1)
        dhyyudy = grad(2)
        dhyyudz = grad(3)
        call grgrad1g_midpoint(hyzu,grad,irg,itg,ipg)
        dhyzudx = grad(1)
        dhyzudy = grad(2)
        dhyzudz = grad(3)
        call grgrad1g_midpoint(hzzu,grad,irg,itg,ipg)
        dhzzudx = grad(1)
        dhzzudy = grad(2)
        dhzzudz = grad(3)
!
        dhyxudx = dhxyudx
        dhyxudy = dhxyudy
        dhyxudz = dhxyudz
        dhzxudx = dhxzudx
        dhzxudy = dhxzudy
        dhzxudz = dhxzudz
        dhzyudx = dhyzudx
        dhzyudy = dhyzudy
        dhzyudz = dhyzudz
!
! --- end
!
        dhdh = dhxxudx*dhxxddx +  dhxxudy*dhxxddy +  dhxxudz*dhxxddz &
        &    + dhxyudx*dhxyddx +  dhxyudy*dhxyddy +  dhxyudz*dhxyddz &
        &    + dhxzudx*dhxzddx +  dhxzudy*dhxzddy +  dhxzudz*dhxzddz &
        &    + dhyxudx*dhyxddx +  dhyxudy*dhyxddy +  dhyxudz*dhyxddz &
        &    + dhyyudx*dhyyddx +  dhyyudy*dhyyddy +  dhyyudz*dhyyddz &
        &    + dhyzudx*dhyzddx +  dhyzudy*dhyzddy +  dhyzudz*dhyzddz &
        &    + dhzxudx*dhzxddx +  dhzxudy*dhzxddy +  dhzxudz*dhzxddz &
        &    + dhzyudx*dhzyddx +  dhzyudy*dhzyddy +  dhzyudz*dhzyddz &
        &    + dhzzudx*dhzzddx +  dhzzudy*dhzzddy +  dhzzudz*dhzzddz
!
        dphdph = dphxxu*dphxxd + dphxyu*dphxyd + dphxzu*dphxzd &
        &      + dphyxu*dphyxd + dphyyu*dphyyd + dphyzu*dphyzd &
        &      + dphzxu*dphzxd + dphzyu*dphzyd + dphzzu*dphzzd
!
        heli(2,1) = heli(1,2)
        heli(3,1) = heli(1,3)
        heli(3,2) = heli(2,3)
        rlieaij(2,1) = rlieaij(1,2)
        rlieaij(3,1) = rlieaij(1,3)
        rlieaij(3,2) = rlieaij(2,3)
          rnlab(2,1) =   rnlab(1,2)
          rnlab(3,1) =   rnlab(1,3)
          rnlab(3,2) =   rnlab(2,3)
         rpddal(2,1) =  rpddal(1,2)
         rpddal(3,1) =  rpddal(1,3)
         rpddal(3,2) =  rpddal(2,3)
         aacabc(2,1) =  aacabc(1,2)
         aacabc(3,1) =  aacabc(1,3)
         aacabc(3,2) =  aacabc(2,3)
!            sab(2,1) =     sab(1,2)
!            sab(3,1) =     sab(1,3)
!            sab(3,2) =     sab(2,3)
          r2dab(2,1) =   r2dab(1,2)
          r2dab(3,1) =   r2dab(1,3)
          r2dab(3,2) =   r2dab(2,3)
        liebaij(2,1) = liebaij(1,2)
        liebaij(3,1) = liebaij(1,3)
        liebaij(3,2) = liebaij(2,3)
!
        trheli = 0.0d0
        traij = 0.0d0
        trlieaij = 0.0d0
        trrnlab = 0.0d0
        trrpddal = 0.0d0
        traa = 0.0d0
!        trsab = 0.0d0
        trd2ph = 0.0d0
        trr2dab = 0.0d0
        trliebaij = 0.0d0
        do ib = 1, 3
          do ia = 1, 3
            trheli    = trheli    + gamu(ia,ib)*   heli(ia,ib)
            traij     = traij     + gamu(ia,ib)*    aij(ia,ib)
            trlieaij  = trlieaij  + gamu(ia,ib)*rlieaij(ia,ib)
            trrnlab   = trrnlab   + gamu(ia,ib)*  rnlab(ia,ib)
            trrpddal  = trrpddal  + gamu(ia,ib)* rpddal(ia,ib)
            traa      = traa      + gamu(ia,ib)* aacabc(ia,ib)
!            trsab    = trsab     + gamu(ia,ib)*    sab(ia,ib)
            trd2ph    = trd2ph    + gamu(ia,ib)* d2phab(ia,ib)
            trr2dab   = trr2dab   + gamu(ia,ib)*  r2dab(ia,ib)
            trliebaij = trliebaij + gamu(ia,ib)*liebaij(ia,ib)
          end do
        end do
!
        do ic = 1, 6
!
          ia = 1 + ic/4 + ic/6
          ib = ic - (ic/4)*2 - ic/6
!
          sgam = san*gamd(ia,ib)
          tfaij     =     aij(ia,ib) - sgam*traij
          tflieaij  = rlieaij(ia,ib) - sgam*trlieaij
          tfrnlab   =   rnlab(ia,ib) - sgam*trrnlab 
          tfrpddal  =  rpddal(ia,ib) - sgam*trrpddal
          tfaa      =  aacabc(ia,ib) - sgam*traa
          tfr2dab   =   r2dab(ia,ib) - sgam*trr2dab 
          tfliebaij = liebaij(ia,ib) - sgam*trliebaij
          tfd2ph    =  d2phab(ia,ib) - sgam*trd2ph   
!
          gamdhdh  = gamd(ia,ib)*dhdh
          ogdphdph = 0.0d0
          if (chope == 'H') ogdphdph = gamd(ia,ib)*ome**2*dphdph
!
          tfliep4aij = 0.0d0
          tfheli     = 0.0d0

          if (chgra == 'w') tfliep4aij = tfaij*ovpadpsi + ps4oal*tflieaij

          if (chgra == 'k') tfliep4aij = tfaij*bepadpsi + ps4oal*tfliebaij


!          if(itg.eq.ntg/2.and.ipg.eq.1.and.ic.eq.1) &
!         & write(6,*) tfliep4aij, tfaij, ovpadpsi, ps4oal, tflieaij
          if (chgra == 'c') tfliep4aij = tfaij*ovpadpsi + ps4oal*tflieaij
          if (chgra == 'C') tfliep4aij = tfaij*ovpadpsi + ps4oal*tflieaij &
             &                             - 0.5d0*ome**2*tfd2ph*cutoff
          if (chgra == 'h') tfheli = heli(ia,ib) - sgam*trheli
          if (chgra == 'H') tfheli = heli(ia,ib) - sgam*trheli
          if (chgra == 'W') then 
                            tfheli = heli(ia,ib) - sgam*trheli
                            tfheli = tfheli*cutoff
                            tfliep4aij = tfaij*ovpadpsi + ps4oal*tflieaij
                            tfliep4aij = tfliep4aij*(1.0d0 - cutoff)
          end if
!

          souten(irg,itg,ipg,ic) = 2.0d0*( tfheli                                      &
                    &            + tfliep4aij + tfr2dab + tfrnlab + tfrpddal           &
                    &            + psigc4*(trk(irg,itg,ipg)*tfaij/3.0d0 - 2.0d0*tfaa)  &
                    &            - roku*gamdhdh + roku*ogdphdph ) 
!
        end do
      end do
    end do
  end do
!
!
  deallocate(dfdx)
  deallocate(dfdy)
  deallocate(dfdz)
  deallocate(fnc0)
  deallocate(grada)
  deallocate(gradp)
  deallocate(gradap)
  deallocate(gradap2)
  deallocate(gradk)
  deallocate(gradb)
!
end subroutine sourceterm_trfreeG_WL_bhex
