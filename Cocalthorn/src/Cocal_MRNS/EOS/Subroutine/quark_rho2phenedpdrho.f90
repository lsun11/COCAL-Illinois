Subroutine quark_rho2phenedpdrho(rho,pre,hh,ene,dpdrho)
  use def_qeos_parameter
  implicit none
  integer :: ii
  real(8),intent(inout) :: rho, pre, hh, ene, dpdrho
  pre=0.0d0
  hh =0.0d0
  ene=0.0d0
  dpdrho=0.0d0
  if (rho.le.0.0d0) rho=0.0d0
  do ii=1,nphase
     pre=abc(ii)*rho**abi(ii)+pre
     hh =abch(ii)*rho**abih(ii)+hh
     ene=abcene(ii)*rho**abiene(ii)+ene
     dpdrho=abi(ii)*abc(ii)*rho**(abi(ii)-1.0d0)+dpdrho
  end do
  do ii=nphase+1,2*nphase
     hh=abch(ii)*rho**abih(ii)+hh
  end do
     ene=ene+rho+rho*eneconst_gcm1
     hh=hh+1.0d0+eneconst_gcm1
END Subroutine quark_rho2phenedpdrho     
