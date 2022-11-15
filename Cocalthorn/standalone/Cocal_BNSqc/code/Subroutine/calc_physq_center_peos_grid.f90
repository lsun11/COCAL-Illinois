subroutine calc_physq_center_peos_grid
  use grid_parameter, only  :  nrf, ntf, npf, ntfxy
  use def_matter, only : emd
  use def_quantities, only : rho_c, pre_c, epsi_c, q_c, &
  &                          rho_max, pre_max, epsi_max, q_max 
  implicit none
  real(8) :: xsol, hh, emdmax
  integer :: irf,itf,ipf, irfmax, itfmax, ipfmax
!
  q_c = emd(0,0,0)
  call peos_q2hprho(q_c, hh, pre_c, rho_c, epsi_c)
!  write(6,'(a6,1p,e23.15,a21,i3,a1,i3,a1,i3,a1)') "q_c  =", q_c,"   at (irf,itf,ipf)=(", &
!  &                                                0, ",", 0, ",", 0,")"
  write(6,'(a6,1p,e23.15,a13)') "q_c  =", q_c,"   at (0,0,0)"
!
  itf=ntfxy;  ipf=0
  q_max = 0.0d0;  itfmax=ntfxy;  ipfmax=0
  do irf=0, nrf
    if ( emd(irf,itf,ipf) > q_max )  then
      irfmax = irf
      q_max = emd(irf,itf,ipf)
    end if
  end do
!
  if (irfmax .ne. 0) then
    call search_emdmax_xaxis(irfmax, emdmax, xsol) 
    q_max = emdmax
    write(6,'(a6,1p,e23.15,a8)') "q_max=", q_max,"   at x=", xsol 
  else
    write(6,'(a6,1p,e23.15,a13)') "q_max=", q_max,"   at (0,0,0)"
  end if
!  write(6,'(a6,1p,e23.15,a21,i3,a1,i3,a1,i3,a1)') "q_max=", q_max,"   at (irf,itf,ipf)=(", &
!    &                                                   irfmax, ",", itfmax, ",", ipfmax,")"
!
  call peos_q2hprho(q_max, hh, pre_max, rho_max, epsi_max)
!
end subroutine calc_physq_center_peos_grid
