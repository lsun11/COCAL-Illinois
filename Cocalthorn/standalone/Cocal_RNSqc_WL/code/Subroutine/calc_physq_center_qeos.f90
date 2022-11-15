subroutine calc_physq_center_qeos
  use grid_parameter, only  :  nrf, ntf, npf
  use def_matter, only : rhof 
  use def_quantities, only : rho_c, pre_c, epsi_c, q_c, &
  &                          rho_max, pre_max, epsi_max, q_max 
  implicit none
  real(8) :: hh, emdmax=-1.0d0, dummy
  integer :: irf,itf,ipf, irfmax, itfmax, ipfmax
!
  rho_c = rhof(0,0,0)
  call quark_rho2phenedpdrho(rho_c, pre_c, hh, epsi_c, dummy)
  q_c = pre_c/rho_c
  write(6,'(a6,1p,e23.15,a21,i3,a1,i3,a1,i3,a1)') "q_c  =", q_c,"   at (irf,itf,ipf)=(", &
  &                                                0, ",", 0, ",", 0,")"
!
  rho_max = 0.0d0
  do irf=0, nrf
    do itf=0, ntf
      do ipf=0, npf
        if ( rhof(irf,itf,ipf) > rho_max )  then
          irfmax = irf
          itfmax = itf
          ipfmax = ipf
          rho_max = rhof(irf,itf,ipf)
        end if
      end do
    end  do 
  end do
!
  call quark_rho2phenedpdrho(rho_max, pre_max, hh, epsi_max, dummy)
  q_max = pre_max/rho_max
  write(6,'(a6,1p,e23.15,a21,i3,a1,i3,a1,i3,a1)') "q_max=", q_max,"   at (irf,itf,ipf)=(", &
  &                                                irfmax, ",", itfmax, ",", ipfmax,")"
!
!
end subroutine calc_physq_center_qeos
