subroutine quark_sound_speed(c_s,rho0)
  use def_qeos_parameter
  implicit none
  integer :: ii
  real(8), intent(inout) :: c_s, rho0
  real(8) :: dpdrho, denedrho
  ii=0
  dpdrho=0.0d0
  denedrho=1.0d0 + eneconst_gcm1
  c_s=0
  do ii=1, nphase
     dpdrho=abi(ii)*abc(ii)*rho0**(abi(ii)-1.0d0)+dpdrho
     denedrho=(abc(ii)*abi(ii)/(abi(ii)-1.0d0))*rho0**(abi(ii)-1.0d0)+denedrho
  end do
  c_s=(dpdrho/denedrho)**(1./2.)
!  if (c_s.gt.1) then
!     write(6,*) 'superluminal sound speed', c_s
!  end if 
  
end subroutine quark_sound_speed
