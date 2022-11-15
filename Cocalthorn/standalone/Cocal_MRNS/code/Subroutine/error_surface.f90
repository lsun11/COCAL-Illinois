subroutine error_surface(pot,pot_bak,error,flag)
  use phys_constant,  only : long
  use grid_parameter, only : nrg, ntg, npg, eps_coc
  implicit none
  real(long), pointer :: pot(:,:), pot_bak(:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: itg, ipg
!
! The criteria for error in rs is smaller than the other variables.
  error = 0.0d0
  flag = 0
  do itg = 0, ntg
    do ipg = 0, npg
      error_pot = 2.0d0*dabs(pot(itg,ipg) -     pot_bak(itg,ipg)) &
    &                 /(dabs(pot(itg,ipg))+dabs(pot_bak(itg,ipg))+small)
      if (error_pot > eps_coc*0.3) flag = 1
      error = dmax1(error,error_pot)
    end do
  end do
end subroutine error_surface
