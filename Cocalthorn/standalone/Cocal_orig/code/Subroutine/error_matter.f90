subroutine error_matter(pot,pot_bak,error,flag)
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntf, npf, eps
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long), intent(out) :: error
  integer,    intent(out) :: flag
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irf, itf, ipf
!
  error = 0.0d0
  flag = 0
  do irf = 0, nrf-1
    do itf = 0, ntf
      do ipf = 0, npf
        error_pot = 2.0d0*dabs(pot(irf,itf,ipf) -     pot_bak(irf,itf,ipf)) &
      &                 /(dabs(pot(irf,itf,ipf))+dabs(pot_bak(irf,itf,ipf)) &
      &                 + small)
        if (error_pot > eps) flag = 1
        error = dmax1(error,error_pot)
      end do
    end do
  end do
end subroutine error_matter
