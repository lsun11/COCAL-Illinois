Subroutine mit_ene2prho(ene,pre,rho)
  use def_qeos_parameter
  implicit none
  real(8),intent(inout) :: rho, pre, ene
  pre=aq*(ene-enesurf_gcm1)
  rho=0.0d0 
END Subroutine mit_ene2prho   
