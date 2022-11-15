subroutine error_adjust_parameter(niq,pot,pot_bak,error,flag)
  use phys_constant,  only : long
  use grid_parameter, only : eps_coc
  implicit none
  real(long), intent(in)  :: pot(niq), pot_bak(niq)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer :: niq, ii
!
  error = 0.0d0
  flag = 0
  do ii = 1, niq
    error_pot = 2.0d0*dabs(pot(ii) -     pot_bak(ii)) &
    &               /(dabs(pot(ii))+dabs(pot_bak(ii)) + small)
    if (error_pot > eps_coc) flag = 1
    error = dmax1(error,error_pot)
  end do
end subroutine error_adjust_parameter
