Function quark_rho2p(rho) 
  use def_qeos_parameter, only: nphase, abc, abi
  implicit none
  real(8) :: rho, quark_rho2p
  integer :: ii
  quark_rho2p = 0.0d0

!  write(6,'(a20,i2,1p,2e23.15)') "abc,rho", nphase, abc(1), rho
 
  do ii=1,nphase
     quark_rho2p=abc(ii)*rho**abi(ii)+quark_rho2p
  end do

!  write(6,'(a20,1p,2e23.15)') "quark_rho2p=", quark_rho2p

  return
END Function quark_rho2p
