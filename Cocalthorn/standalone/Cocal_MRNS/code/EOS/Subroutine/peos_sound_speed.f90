subroutine peos_sound_speed(q,ss)
!
  use def_peos_parameter	!abc,abi,rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: ss
  real(8)                :: hh, pre, rho, epsi, abin
  integer                :: iphase
!
  call peos_q2hprho(q, hh, pre, rho, epsi)
  call peos_lookup(q, qi, iphase)
  abin = abi(iphase)
  ss = sqrt(abin*pre/rho/hh)
!
end subroutine peos_sound_speed
