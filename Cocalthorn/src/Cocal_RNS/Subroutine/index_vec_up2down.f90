subroutine index_vec_up2down(vecxd,vecyd,veczd,vecxu,vecyu,veczu)
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_hij, only : hxxd, hxyd, hxzd, hyyd, hyzd, hzzd
  implicit none
  real(8), pointer :: vecxu(:,:,:), vecyu(:,:,:), veczu(:,:,:), &
  &                   vecxd(:,:,:), vecyd(:,:,:), veczd(:,:,:)
  real(8) :: gmxxd, gmxyd, gmxzd, gmyxd, gmyyd, gmyzd, &
  &          gmzxd, gmzyd, gmzzd
  integer :: ipg, itg, irg
!
! --- Lowering index of shift.  
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
end subroutine index_vec_up2down
