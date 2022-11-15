subroutine peos_h2qprho(h,q,pre,rho,ened)
!
  use def_peos_parameter	!abc,abi,rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(in)  :: h
  real(8), intent(out) :: q, pre, rho, ened
  real(8)              :: qq, hin, qin, abin, abct
  real(8)              :: fac1, fac2, fack, small = 1.0d-16
  integer              :: iphase
!
  call peos_lookup(h, hi, iphase)
  hin  = hi(iphase)
  qin  = qi(iphase)
  abin = abi(iphase)
  abct = abc(iphase)
!
  fac1 = 1.0d0/(abin - 1.0d0)
  fac2 = abin/(abin - 1.0d0)
  fack = abct**(-fac1)
!
  q = (h - hin)/fac2 + qin
  qq = q
  if (h <= 1.0d0) qq = small/fac2
  pre = fack*qq**fac2
  rho = fack*qq**fac1
  ened = rho*h - pre
!
end subroutine peos_h2qprho
