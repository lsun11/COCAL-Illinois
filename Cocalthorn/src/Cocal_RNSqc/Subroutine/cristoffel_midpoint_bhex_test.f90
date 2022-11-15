subroutine cristoffel_midpoint_bhex
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd, &
  &                          hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_metric
  use def_cristoffel, only : cri
  use def_gamma_crist, only : gmcrix, gmcriy, gmcriz
  use def_formulation, only : chgra
  use def_kerr_schild, only: Fxu_ks, Fyu_ks, Fzu_ks
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use make_array_3d
  implicit none
  real(8) :: crid11, crid12, crid13, crid14, crid15, crid16, &
  &          crid21, crid22, crid23, crid24, crid25, crid26, &
  &          crid31, crid32, crid33, crid34, crid35, crid36, &
  &          dhyxdx, dhyxdy, dhyxdz, dhzxdx, dhzxdy, dhzxdz, &
  &          dhzydx, dhzydy, dhzydz, &
  &          dhxxdx, dhxxdy, dhxxdz, dhxydx, dhxydy, dhxydz, &
  &          dhxzdx, dhxzdy, dhxzdz, dhyydx, dhyydy, dhyydz, &
  &          dhyzdx, dhyzdy, dhyzdz, dhzzdx, dhzzdy, dhzzdz, &
  &          gmxxu, gmxyu, gmxzu, gmyyu, gmyzu, gmzzu, &
  &          gmyxu, gmzxu, gmzyu, ps4oal2, bbxx, bbxy, bbxz, bbyy, &
  &          bbyz, bbzz, psim, alphm, bvxum, bvyum, bvzum

  real(8), pointer :: dfxxdx(:,:,:), dfxxdy(:,:,:), dfxxdz(:,:,:), &
  &                   dfxydx(:,:,:), dfxydy(:,:,:), dfxydz(:,:,:), &
  &                   dfxzdx(:,:,:), dfxzdy(:,:,:), dfxzdz(:,:,:), &
  &                   dfyydx(:,:,:), dfyydy(:,:,:), dfyydz(:,:,:), &
  &                   dfyzdx(:,:,:), dfyzdy(:,:,:), dfyzdz(:,:,:), &
  &                   dfzzdx(:,:,:), dfzzdy(:,:,:), dfzzdz(:,:,:)
  integer :: ipg, itg, irg
!
  call alloc_array3d(dfxxdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxxdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxxdz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxydx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxydy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxydz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxzdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxzdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfxzdz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfyydx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfyydy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfyydz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfyzdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfyzdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfyzdz,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfzzdx,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfzzdy,1,nrg,1,ntg,1,npg)
  call alloc_array3d(dfzzdz,1,nrg,1,ntg,1,npg)
!
! --- Compute Cristoffel symbols,   
!     whose value is assigned on the grid points. 
!     C^a_bc=1/2 gm^ad(pa_b h_dc+...)
! --- for interporation.
!
  call grgrad_midpoint(hxxd,dfxxdx,dfxxdy,dfxxdz)
  call grgrad_midpoint(hxyd,dfxydx,dfxydy,dfxydz)
  call grgrad_midpoint(hxzd,dfxzdx,dfxzdy,dfxzdz)
  call grgrad_midpoint(hyyd,dfyydx,dfyydy,dfyydz)
  call grgrad_midpoint(hyzd,dfyzdx,dfyzdy,dfyzdz)
  call grgrad_midpoint(hzzd,dfzzdx,dfzzdy,dfzzdz)
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
!
        call interpo_linear_type0(psim, psi, irg,itg,ipg)
        call interpo_linear_type0(alphm,alph,irg,itg,ipg)
        call interpo_linear_type0(bvxum,bvxu,irg,itg,ipg)
        call interpo_linear_type0(bvyum,bvyu,irg,itg,ipg)
        call interpo_linear_type0(bvzum,bvzu,irg,itg,ipg)

        call interpo_linear_type0(gmxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(gmxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(gmxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(gmyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(gmyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(gmzzu,hzzu,irg,itg,ipg)

!        ps4oal2 = psim**4/alphm**2
!        bbxx = ps4oal2*bvxum*bvxum
!        bbxy = ps4oal2*bvxum*bvyum
!        bbxz = ps4oal2*bvxum*bvzum
!        bbyy = ps4oal2*bvyum*bvyum
!        bbyz = ps4oal2*bvyum*bvzum
!        bbzz = ps4oal2*bvzum*bvzum

        gmxxu = gmxxu + 1.0d0 !+ bbxx
        gmyyu = gmyyu + 1.0d0 !+ bbyy
        gmzzu = gmzzu + 1.0d0 !+ bbzz
        gmxyu = gmxyu         !+ bbxy
        gmxzu = gmxzu         !+ bbxz
        gmyzu = gmyzu         !+ bbyz
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu
!
        dhxxdx = dfxxdx(irg,itg,ipg)
        dhxxdy = dfxxdy(irg,itg,ipg)
        dhxxdz = dfxxdz(irg,itg,ipg)
!
        dhxydx = dfxydx(irg,itg,ipg)
        dhxydy = dfxydy(irg,itg,ipg)
        dhxydz = dfxydz(irg,itg,ipg)
!
        dhxzdx = dfxzdx(irg,itg,ipg)
        dhxzdy = dfxzdy(irg,itg,ipg)
        dhxzdz = dfxzdz(irg,itg,ipg)
!
        dhyydx = dfyydx(irg,itg,ipg)
        dhyydy = dfyydy(irg,itg,ipg)
        dhyydz = dfyydz(irg,itg,ipg)
!
        dhyzdx = dfyzdx(irg,itg,ipg)
        dhyzdy = dfyzdy(irg,itg,ipg)
        dhyzdz = dfyzdz(irg,itg,ipg)
!
        dhzzdx = dfzzdx(irg,itg,ipg)
        dhzzdy = dfzzdy(irg,itg,ipg)
        dhzzdz = dfzzdz(irg,itg,ipg)
!
        dhyxdx = dhxydx
        dhyxdy = dhxydy
        dhyxdz = dhxydz
!
        dhzxdx = dhxzdx
        dhzxdy = dhxzdy
        dhzxdz = dhxzdz
!
        dhzydx = dhyzdx
        dhzydy = dhyzdy
        dhzydz = dhyzdz
!
! --- Impose Dirac gauge.
!
!dirac      dhx = gmxxu*dhxxdx + gmxyu*dhxxdy + gmxzu*dhxxdz
!dirac     &    + gmyxu*dhyxdx + gmyyu*dhyxdy + gmyzu*dhyxdz
!dirac     &    + gmzxu*dhzxdx + gmzyu*dhzxdy + gmzzu*dhzxdz
!dirac      dhy = gmxxu*dhxydx + gmxyu*dhxydy + gmxzu*dhxydz
!dirac     &    + gmyxu*dhyydx + gmyyu*dhyydy + gmyzu*dhyydz
!dirac     &    + gmzxu*dhzydx + gmzyu*dhzydy + gmzzu*dhzydz
!dirac      dhz = gmxxu*dhxzdx + gmxyu*dhxzdy + gmxzu*dhxzdz
!dirac     &    + gmyxu*dhyzdx + gmyyu*dhyzdy + gmyzu*dhyzdz
!dirac     &    + gmzxu*dhzzdx + gmzyu*dhzzdy + gmzzu*dhzzdz
!dirac      gzzinv = 1.0d0/gmzzu
!dirac      dhxzdz = - gzzinv*dhx + dhxzdz
!dirac      dhyzdz = - gzzinv*dhy + dhyzdz
!dirac      dhzzdz = - gzzinv*dhz + dhzzdz
!dirac      dhzxdz = dhxzdz
!dirac      dhzydz = dhyzdz
!
! --- end
!
        crid11 = 0.5d0*(dhxxdx + dhxxdx - dhxxdx)
        crid12 = 0.5d0*(dhxxdy + dhxydx - dhxydx)
        crid13 = 0.5d0*(dhxxdz + dhxzdx - dhxzdx)
        crid14 = 0.5d0*(dhxydy + dhxydy - dhyydx)
        crid15 = 0.5d0*(dhxydz + dhxzdy - dhyzdx)
        crid16 = 0.5d0*(dhxzdz + dhxzdz - dhzzdx)
        crid21 = 0.5d0*(dhxydx + dhxydx - dhxxdy)
        crid22 = 0.5d0*(dhxydy + dhyydx - dhxydy)
        crid23 = 0.5d0*(dhxydz + dhyzdx - dhxzdy)
        crid24 = 0.5d0*(dhyydy + dhyydy - dhyydy)
        crid25 = 0.5d0*(dhyydz + dhyzdy - dhyzdy)
        crid26 = 0.5d0*(dhyzdz + dhyzdz - dhzzdy)
        crid31 = 0.5d0*(dhxzdx + dhxzdx - dhxxdz)
        crid32 = 0.5d0*(dhxzdy + dhyzdx - dhxydz)
        crid33 = 0.5d0*(dhxzdz + dhzzdx - dhxzdz)
        crid34 = 0.5d0*(dhyzdy + dhyzdy - dhyydz)
        crid35 = 0.5d0*(dhyzdz + dhzzdy - dhyzdz)
        crid36 = 0.5d0*(dhzzdz + dhzzdz - dhzzdz)
!
        cri(irg,itg,ipg,1,1) = &
           &      (-gmxzu*dhxxdz - gmxyu*dhxxdy + gmxxu*dhxxdx + &
           &  2.0d0*gmxyu*dhxydx + 2.0d0*gmxzu*dhxzdx)*0.5d0
        cri(irg,itg,ipg,1,2) = &
           &      (-gmxzu*dhxydz + gmxxu*dhxxdy + gmxzu*dhxzdy + &
           &        gmxyu*dhyydx + gmxzu*dhyzdx)*0.5d0
        cri(irg,itg,ipg,1,3) = &
           &       (gmxxu*dhxxdz + gmxyu*dhxydz - gmxyu*dhxzdy + &
           &        gmxyu*dhyzdx + gmxzu*dhzzdx)*0.5d0
        cri(irg,itg,ipg,1,4) = &
           &      (-gmxzu*dhyydz + 2.0d0*gmxxu*dhxydy + gmxyu*dhyydy + &
           &  2.0d0*gmxzu*dhyzdy - gmxxu*dhyydx)*0.5d0
        cri(irg,itg,ipg,1,5) = &
           &       (gmxxu*dhxydz + gmxyu*dhyydz + gmxxu*dhxzdy + &
           &        gmxzu*dhzzdy - gmxxu*dhyzdx)*0.5d0
        cri(irg,itg,ipg,1,6) = &
           & (2.0d0*gmxxu*dhxzdz + 2.0d0*gmxyu*dhyzdz + &
           &        gmxzu*dhzzdz - gmxyu*dhzzdy - gmxxu*dhzzdx)*0.5d0
!
        cri(irg,itg,ipg,2,1) = &
           &      (-gmyzu*dhxxdz - gmyyu*dhxxdy + gmxyu*dhxxdx + &
           &  2.0d0*gmyyu*dhxydx + 2.0d0*gmyzu*dhxzdx)*0.5d0
        cri(irg,itg,ipg,2,2) = &
           &      (-gmyzu*dhxydz + gmxyu*dhxxdy + gmyzu*dhxzdy + &
           &        gmyyu*dhyydx + gmyzu*dhyzdx)*0.5d0
        cri(irg,itg,ipg,2,3) = &
           &       (gmxyu*dhxxdz + gmyyu*dhxydz - gmyyu*dhxzdy + &
           &        gmyyu*dhyzdx + gmyzu*dhzzdx)*0.5d0
        cri(irg,itg,ipg,2,4) = &
           &      (-gmyzu*dhyydz + 2.0d0*gmxyu*dhxydy + gmyyu*dhyydy + &
           &  2.0d0*gmyzu*dhyzdy - gmxyu*dhyydx)*0.5d0
        cri(irg,itg,ipg,2,5) = &
           &       (gmxyu*dhxydz + gmyyu*dhyydz + gmxyu*dhxzdy + &
           &        gmyzu*dhzzdy - gmxyu*dhyzdx)*0.5d0
        cri(irg,itg,ipg,2,6) = &
           & (2.0d0*gmxyu*dhxzdz + 2.0d0*gmyyu*dhyzdz + &
           &        gmyzu*dhzzdz - gmyyu*dhzzdy - gmxyu*dhzzdx)*0.5d0
!
        cri(irg,itg,ipg,3,1) = &
           &      (-gmzzu*dhxxdz - gmyzu*dhxxdy + gmxzu*dhxxdx + &
           &  2.0d0*gmyzu*dhxydx + 2.0d0*gmzzu*dhxzdx)*0.5d0
        cri(irg,itg,ipg,3,2) = &
           &      (-gmzzu*dhxydz + gmxzu*dhxxdy + gmzzu*dhxzdy + &
           &        gmyzu*dhyydx + gmzzu*dhyzdx)*0.5d0
        cri(irg,itg,ipg,3,3) = &
           &       (gmxzu*dhxxdz + gmyzu*dhxydz - gmyzu*dhxzdy + &
           &        gmyzu*dhyzdx + gmzzu*dhzzdx)*0.5d0
        cri(irg,itg,ipg,3,4) = &
           &      (-gmzzu*dhyydz + 2.0d0*gmxzu*dhxydy + gmyzu*dhyydy + &
           &  2.0d0*gmzzu*dhyzdy - gmxzu*dhyydx)*0.5d0
        cri(irg,itg,ipg,3,5) = &
           &       (gmxzu*dhxydz + gmyzu*dhyydz + gmxzu*dhxzdy + &
           &        gmzzu*dhzzdy - gmxzu*dhyzdx)*0.5d0
        cri(irg,itg,ipg,3,6) = &
           & (2.0d0*gmxzu*dhxzdz + 2.0d0*gmyzu*dhyzdz + &
           &        gmzzu*dhzzdz - gmyzu*dhzzdy - gmxzu*dhzzdx)*0.5d0
!
! --  gamma^ab C^c_ab
!
!gauge      gmcrix(irg,itg,ipg) = gmxxu*cri(irg,itg,ipg,1,1)
!gauge     &                    + gmxyu*cri(irg,itg,ipg,1,2)*2.0d0
!gauge     &                    + gmxzu*cri(irg,itg,ipg,1,3)*2.0d0
!gauge     &                    + gmyyu*cri(irg,itg,ipg,1,4)
!gauge     &                    + gmyzu*cri(irg,itg,ipg,1,5)*2.0d0
!gauge     &                    + gmzzu*cri(irg,itg,ipg,1,6)
!gauge      gmcriy(irg,itg,ipg) = gmxxu*cri(irg,itg,ipg,2,1)
!gauge     &                    + gmxyu*cri(irg,itg,ipg,2,2)*2.0d0
!gauge     &                    + gmxzu*cri(irg,itg,ipg,2,3)*2.0d0
!gauge     &                    + gmyyu*cri(irg,itg,ipg,2,4)
!gauge     &                    + gmyzu*cri(irg,itg,ipg,2,5)*2.0d0
!gauge     &                    + gmzzu*cri(irg,itg,ipg,2,6)
!gauge      gmcriz(irg,itg,ipg) = gmxxu*cri(irg,itg,ipg,3,1)
!gauge     &                    + gmxyu*cri(irg,itg,ipg,3,2)*2.0d0
!gauge     &                    + gmxzu*cri(irg,itg,ipg,3,3)*2.0d0
!gauge     &                    + gmyyu*cri(irg,itg,ipg,3,4)
!gauge     &                    + gmyzu*cri(irg,itg,ipg,3,5)*2.0d0
!gauge     &                    + gmzzu*cri(irg,itg,ipg,3,6)
!
        gmcrix(irg,itg,ipg) = 0.0d0
        gmcriy(irg,itg,ipg) = 0.0d0
        gmcriz(irg,itg,ipg) = 0.0d0
!
        if (chgra == 'k') then
          gmcrix(irg,itg,ipg) = -Fxu_ks(irg,itg,ipg)
          gmcriy(irg,itg,ipg) = -Fyu_ks(irg,itg,ipg)
          gmcriz(irg,itg,ipg) = -Fzu_ks(irg,itg,ipg)
        end if
!
      end do
    end do
  end do
!
!gauge      write(6,*) 'gmcri',gmcrix(8,7,12), gmcriy(8,7,12), gmcriz(8,7,12)
  deallocate(dfxxdx)
  deallocate(dfxxdy)
  deallocate(dfxxdz)
  deallocate(dfxydx)
  deallocate(dfxydy)
  deallocate(dfxydz)
  deallocate(dfxzdx)
  deallocate(dfxzdy)
  deallocate(dfxzdz)
  deallocate(dfyydx)
  deallocate(dfyydy)
  deallocate(dfyydz)
  deallocate(dfyzdx)
  deallocate(dfyzdy)
  deallocate(dfyzdz)
  deallocate(dfzzdx)
  deallocate(dfzzdy)
  deallocate(dfzzdz)
!
end subroutine cristoffel_midpoint_bhex
