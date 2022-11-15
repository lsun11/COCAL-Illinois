subroutine index_vec_down2up_midpoint(vecxu,vecyu,veczu,vecxd,vecyd,veczd)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxu, hxyu, hxzu, hyyu, hyzu, hzzu
  use interface_interpo_linear_type0
  implicit none
  real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
  &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
  real(8) :: hhxxu, hhxyu, hhxzu, hhyyu, hhyzu, hhzzu
  real(8) :: gmxxu, gmxyu, gmxzu, gmyxu, gmyyu, gmyzu, &
  &          gmzxu, gmzyu, gmzzu
  integer :: ipg, itg, irg
!
! --- Rasing index of vector at midpoint.  
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(hhxxu,hxxu,irg,itg,ipg)
        call interpo_linear_type0(hhxyu,hxyu,irg,itg,ipg)
        call interpo_linear_type0(hhxzu,hxzu,irg,itg,ipg)
        call interpo_linear_type0(hhyyu,hyyu,irg,itg,ipg)
        call interpo_linear_type0(hhyzu,hyzu,irg,itg,ipg)
        call interpo_linear_type0(hhzzu,hzzu,irg,itg,ipg)
        gmxxu = 1.0d0 + hhxxu
        gmxyu =         hhxyu
        gmxzu =         hhxzu
        gmyxu =         hhxyu
        gmyyu = 1.0d0 + hhyyu
        gmyzu =         hhyzu
        gmzxu =         hhxzu
        gmzyu =         hhyzu
        gmzzu = 1.0d0 + hhzzu
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
end subroutine index_vec_down2up_midpoint
!bcd2u
