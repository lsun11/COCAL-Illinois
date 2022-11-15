subroutine index_vec_up2down_midpoint(vecxd,vecyd,veczd,vecxu,vecyu,veczu)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  use interface_interpo_linear_type0
  implicit none
  real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
  &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
  real(8) :: hhxxd, hhxyd, hhxzd, hhyyd, hhyzd, hhzzd
  real(8) :: gmxxd, gmxyd, gmxzd, gmyxd, gmyyd, gmyzd, &
  &          gmzxd, gmzyd, gmzzd
  integer :: ipg, itg, irg
!
! --- Rasing index of vector at midpoint.  
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(hhxxd,hxxd,irg,itg,ipg)
        call interpo_linear_type0(hhxyd,hxyd,irg,itg,ipg)
        call interpo_linear_type0(hhxzd,hxzd,irg,itg,ipg)
        call interpo_linear_type0(hhyyd,hyyd,irg,itg,ipg)
        call interpo_linear_type0(hhyzd,hyzd,irg,itg,ipg)
        call interpo_linear_type0(hhzzd,hzzd,irg,itg,ipg)
        gmxxd = 1.0d0 + hhxxd
        gmxyd =         hhxyd
        gmxzd =         hhxzd
        gmyxd =         hhxyd
        gmyyd = 1.0d0 + hhyyd
        gmyzd =         hhyzd
        gmzxd =         hhxzd
        gmzyd =         hhyzd
        gmzzd = 1.0d0 + hhzzd
        vecxd(irg,itg,ipg) = gmxxd*vecxu(irg,itg,ipg) &
        &                  + gmxyd*vecyu(irg,itg,ipg) &
        &                  + gmxzd*veczu(irg,itg,ipg)
        vecyd(irg,itg,ipg) = gmyxd*vecxu(irg,itg,ipg) &
        &                  + gmyyd*vecyu(irg,itg,ipg) &
        &                  + gmyzd*veczu(irg,itg,ipg)
        veczd(irg,itg,ipg) = gmzxd*vecxu(irg,itg,ipg) &
        &                  + gmzyd*vecyu(irg,itg,ipg) &
        &                  + gmzzd*veczu(irg,itg,ipg)
      end do
    end do
  end do
!
end subroutine index_vec_up2down_midpoint
