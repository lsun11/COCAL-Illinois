subroutine compute_wme2psi(pot)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg
  use def_metric_pBH, only : index_wme
  implicit none
  real(long), pointer    :: pot(:,:,:)
  real(long) :: index
  integer :: irg, itg, ipg
!
  index = - 1.0d0/dble(index_wme)
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
         if(pot(irg,itg,ipg).le.0.0d0) then
            pot(irg,itg,ipg) =  1.0d-14
         end if
         pot(irg,itg,ipg) = pot(irg, itg, ipg)**index
      end do
    end do
  end do
!
end subroutine compute_wme2psi
