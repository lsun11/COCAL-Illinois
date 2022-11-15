subroutine peos_q2hprho(q,h,pre,rho,ened)
!
  use def_peos_parameter	!abc,abi,rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: h, pre, rho, ened
  real(8)                :: hin, qin, abin, abct, fac1, fac2, fack, small
  real(8)                :: rhob, frho, dfrho, error
  integer                :: iphase
!
  call peos_lookup(q, qi, iphase)

  if (iphase==(nphase+1)) then
!s1    rho  = (q+dsqrt(q*q-cbar*constqc))/cbar    ! only for sgma = 1.0
!
!
!   For general sigma.
!sg    rho = rhoi(nphase+1)
!sg    do
!sg      rhob  = rho
!sg      frho  = sgma*cbar*rho**(sgma+1.0d0) - (sgma+1.0d0)*q*rho + constqc
!sg      dfrho = sgma*(sgma+1.0d0)*cbar*rho**sgma - (sgma+1.0d0)*q
!sg      rho   = rho - frho/dfrho
!sg      error = 2.d0*(rho - rhob)/(rho + rhob)
!sg      if (abs(error) <= 1.0d-14) exit
!sg    end do
!
!
    pre  = (sgma*cbar*rho**(sgma+1.0d0) + constqc)/(sgma+1.0d0)
    h    = cbar*rho**sgma
    ened = rho*h - pre
  else
    hin  = hi(iphase)
    qin  = qi(iphase)
    abin = abi(iphase)
    abct = abc(iphase)
!
    fac1 = 1.0d0/(abin - 1.0d0)
    fac2 = abin/(abin - 1.0d0)
    fack = abct**(-fac1)
!
    small = 1.0d-60
    if (q <= small) q = small
    h = hin + fac2*(q - qin)
    if (h <= 1.0d0) h = 1.0d0
    pre = fack*q**fac2
    rho = fack*q**fac1
    ened = rho*h - pre
  end if
!
end subroutine peos_q2hprho
