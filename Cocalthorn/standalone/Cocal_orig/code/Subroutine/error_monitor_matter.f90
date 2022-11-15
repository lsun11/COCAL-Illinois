subroutine error_monitor_matter(pot,pot_bak,char,irmo,itmo,ipmo)
  use phys_constant,  only : long
  use grid_parameter, only : nrf, ntf, npf, eps
  implicit none
  real(long), pointer :: pot(:,:,:), pot_bak(:,:,:)
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irf, itf, ipf, irmo, itmo, ipmo
  character(len=5) :: char
!
  do irf = 0, nrf-1
    do itf = 0, ntf
      do ipf = 0, npf
        error_pot = 2.0d0*dabs(pot(irf,itf,ipf) -     pot_bak(irf,itf,ipf)) &
      &                 /(dabs(pot(irf,itf,ipf))+dabs(pot_bak(irf,itf,ipf)) &
      &                 + small)
        if(irf.eq.irmo.and.itf.eq.itmo.and.ipf.eq.ipmo) then
          write(6,'(a13,1x,a5,1p,e14.6,3i5)')'error monitor', char, &
          &                                   error_pot, irmo, itmo, ipmo
        end if
      end do
    end do
  end do
end subroutine error_monitor_matter
