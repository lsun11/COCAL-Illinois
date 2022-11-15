subroutine error_monitor_surface(pot,pot_bak,char,itmo,ipmo)
  use phys_constant,  only : long
  use grid_parameter, only : ntf, npf, eps
  implicit none
  real(long), pointer :: pot(:,:), pot_bak(:,:)
  real(long) :: error_pot = 0.0d0, small = 1.0d-14
  integer    :: irf, itf, ipf, irmo, itmo, ipmo
  character(len=5) :: char
!
  do itf = 0, ntf
    do ipf = 0, npf
      error_pot = 2.0d0*dabs(pot(itf,ipf) -     pot_bak(itf,ipf)) &
    &                 /(dabs(pot(itf,ipf))+dabs(pot_bak(itf,ipf)) &
    &                 + small)
      if(itf.eq.itmo.and.ipf.eq.ipmo) then
        write(6,'(a13,1x,a5,1p,e14.6,2i5)')'error monitor', char, &
        &                                   error_pot, itmo, ipmo
      end if
    end do
  end do
end subroutine error_monitor_surface
