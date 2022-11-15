subroutine hydro_irbns_1stint_WL_peos(emd)
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC, only : alphf, psif &
  &                             hxxuf, hxyuf, hxzuf, hyyuf, hyzuf, hzzuf
  use make_array_3d
  use interface_interpo_gr2fl
  implicit none
  real(long), pointer :: emd(:,:,:)
  real(long) :: hh, ut, pre, rho, ene, qq
  real(long) :: psic, alphc, lamc
  real(long) :: gamxxu, gamxyu, gamxzu, gamyxu, gamyyu, gamyzu, &
  &             gamzxu, gamzyu, gamzzu
  real(long) :: dxvp, dyvp, dzvp
  integer    :: ir, it, ip
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
!
        gamxxu = hxxuf(ir,it,ip) + 1.0d0
        gamxyu = hxyuf(ir,it,ip)
        gamxzu = hxzuf(ir,it,ip)
        gamyyu = hyyuf(ir,it,ip) + 1.0d0
        gamyzu = hyzuf(ir,it,ip)
        gamzzu = hzzuf(ir,it,ip) + 1.0d0
        gamyxu = gamxyu
        gamzxu = gamxzu
        gamzyu = gamyzu
        dxvp = grad_velpot(ir,it,ip,1)
        dyvp = grad_velpot(ir,it,ip,2)
        dzvp = grad_velpot(ir,it,ip,3)
!
        psic  = psif(ir,it,ip)
        alphc = alphf(ir,it,ip)
        lamc  = lambda(ir,it,ip)
!
        hh = &
       &  ((lamc/alphc)**2 - 1.0d0/psic**4* &
       &  (gamxxu*dxvp*dxvp + gamxyu*dxvp*dyvp + gamxzu*dxvp*dzvp &
       & + gamyxu*dyvp*dxvp + gamyyu*dyvp*dyvp + gamyzu*dyvp*dzvp &
       & + gamzxu*dzvp*dxvp + gamzyu*dzvp*dyvp + gamzzu*dzvp*dzvp))**0.5
!
        call peos_h2qprho(hh, qq, pre, rho, ene)
        emd(ir,it,ip) = qq
!
      end do
    end do
  end do
!
end subroutine hydro_irbns_1stint_WL_peos
