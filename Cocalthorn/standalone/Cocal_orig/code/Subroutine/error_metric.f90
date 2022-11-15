subroutine error_metric(pot,pot_bak,error,flag)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, eps
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irg, itg, ipg
!
  error = 0.0d0
  flag = 0
  do irg = 0, nrg-1
    do itg = 0, ntg
      do ipg = 0, npg
        error_pot = 2.0d0*dabs(pot(irg,itg,ipg) -     pot_bak(irg,itg,ipg)) &
      &                 /(dabs(pot(irg,itg,ipg))+dabs(pot_bak(irg,itg,ipg)) &
      &                 + small)
        if (error_pot > eps) flag = 1
        error = dmax1(error,error_pot)
      end do
    end do
  end do
end subroutine error_metric
