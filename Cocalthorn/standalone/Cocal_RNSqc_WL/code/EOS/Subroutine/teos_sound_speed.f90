subroutine peos_sound_speed(q,ss)
!
  use def_peos_parameter	!abc,abi,rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: ss
  real(8)                :: h1, pre1, rho1, ene1
  integer                :: iphase, i0
  real(8)                :: x3(3), fx3(3)
  real(8), external      :: dfdx_2nd

!
  call peos_q2hprho(q, h1, pre1, rho1, ene1)
!  call peos_lookup(q, qi, iphase)
!  abin = abi(iphase)
!  ss = sqrt(abin*pre/rho/hh)

  call peos_lookup(ene1, enei, iphase)
  i0 = min0(max0(iphase-1,0),nphase-2)
  x3(1:3) = enei(i0:i0+2)
  fx3(1:3) = prei(i0:i0+2)

  ss = sqrt(dfdx_2nd(x3,fx3,ene1))
!
end subroutine peos_sound_speed
