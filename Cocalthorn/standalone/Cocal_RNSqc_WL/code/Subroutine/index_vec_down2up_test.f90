subroutine index_vec_down2up(vecxu,vecyu,veczu,vecxd,vecyd,veczd)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use def_metric
  implicit none
  real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
  &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
  real(8) :: gmxxu, gmxyu, gmxzu, gmyxu, gmyyu, gmyzu, &
  &          gmzxu, gmzyu, gmzzu, ps4oal2, bbxx, bbxy, bbxz, bbyy, &
  &          bbyz, bbzz

  integer :: ipg, itg, irg
!
! --- Lowering index of shift.  
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        ps4oal2 = psi(irg,itg,ipg)**4/alph(irg,itg,ipg)**2
        bbxx = ps4oal2*bvxu(irg,itg,ipg)*bvxu(irg,itg,ipg)
        bbxy = ps4oal2*bvxu(irg,itg,ipg)*bvyu(irg,itg,ipg)
        bbxz = ps4oal2*bvxu(irg,itg,ipg)*bvzu(irg,itg,ipg)
        bbyy = ps4oal2*bvyu(irg,itg,ipg)*bvyu(irg,itg,ipg)
        bbyz = ps4oal2*bvyu(irg,itg,ipg)*bvzu(irg,itg,ipg)
        bbzz = ps4oal2*bvzu(irg,itg,ipg)*bvzu(irg,itg,ipg)

        gmxxu = 1.0d0 + hxxu(irg,itg,ipg) + bbxx
        gmxyu =         hxyu(irg,itg,ipg) + bbxy
        gmxzu =         hxzu(irg,itg,ipg) + bbxz
        gmyyu = 1.0d0 + hyyu(irg,itg,ipg) + bbyy
        gmyzu =         hyzu(irg,itg,ipg) + bbyz
        gmzzu = 1.0d0 + hzzu(irg,itg,ipg) + bbzz
        gmyxu = gmxyu
        gmzxu = gmxzu
        gmzyu = gmyzu

        vecxu(irg,itg,ipg) = gmxxu*vecxd(irg,itg,ipg) &
        &                  + gmxyu*vecyd(irg,itg,ipg) &
        &                  + gmxzu*veczd(irg,itg,ipg)

        vecyu(irg,itg,ipg) = gmyxu*vecxd(irg,itg,ipg) &
        &                  + gmyyu*vecyd(irg,itg,ipg) &
        &                  + gmyzu*veczd(irg,itg,ipg)

        veczu(irg,itg,ipg) = gmzxu*vecxd(irg,itg,ipg) &
        &                  + gmzyu*vecyd(irg,itg,ipg) &
        &                  + gmzzu*veczd(irg,itg,ipg)
      end do
    end do
  end do
!
!
end subroutine index_vec_down2up
!bcd2u
