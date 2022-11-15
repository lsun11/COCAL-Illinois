subroutine peos_q2hprho(q,h,pre,rho,ened)
  use def_peos_parameter     !rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(inout) :: q
  real(8), intent(out)   :: h, pre, rho, ened
  integer, save          :: iphase
  integer                :: i0, i, ii
  real(8), external      :: lagint_4th, lagint_2nd
  real(8)                :: x4(4), f4(4)
  real(8)                :: hin, qin, abin, abct, fac1, fac2, fack, small
!  real(8) :: small = 1.0d-30 
!
!  if (q < small) then
!    q=small;  h=1.0d0;  rho=small**2.8d0;  pre=q*rho;  ened=rho*h-pre
    
!  if (q < small) then
!    q=small;  h=1.0d0;  rho=small;  pre=q*rho;  ened=rho*h-pre
!    q=0.0d0;  h=1.0d0;  rho=0.0d0;  pre=0.0d0;  ened=0.0d0
!    return
!  end if
  if(q<=qi(0))  then 
    iphase=1

    hin  = hi(0)
    qin  = qi(0)
    abin = gamma_crust
    abct = kappa_crust
    fac1 = 1.0d0/(abin - 1.0d0)
    fac2 = abin/(abin - 1.0d0)
    fack = abct**(-fac1)
    small = 1.0d-60
    if (q <= small) q = small
    h = hin + fac2*(q - qin)
    if (h <= 1.0d0) h = 1.0d0
    pre = fack*q**fac2
    rho = fack*q**fac1
    ened = rho*h - pre
    return
  else
    ii=iphase
    if(q.lt.qi(iphase)) then
      do i = ii, 1, -1
        if( (q.ge.qi(i-1)) .and. (q.le.qi(i)))  then
          iphase = i
          exit
        end if
      end do
    else
      do i = ii, nphase
        if( (q.ge.qi(i-1)) .and. (q.le.qi(i)))  then
          iphase = i
          exit
        end if
      end do
    endif
  end if

!  do i = 1, nphase
!    if( (q.ge.qi(i-1)) .and. (q.le.qi(i)))  then
!      iphase = i
!      exit
!    end if
!  end do
!  write(6,*) "iphase=", iphase
!
  i0 = min0(max0(iphase-2,0),nphase-3)
 
  x4(1:4) = qi(i0:i0+3)
  f4(1:4) = hi(i0:i0+3)
  h       = lagint_4th(x4,f4,q)
  if (h <= 1.0d0) h = 1.0d0

  f4(1:4) = prei(i0:i0+3)
  pre     = lagint_4th(x4,f4,q)

  f4(1:4) = rhoi(i0:i0+3)
  rho     = lagint_4th(x4,f4,q)

  f4(1:4) = enei(i0:i0+3)
  ened    = lagint_4th(x4,f4,q)
!
end subroutine peos_q2hprho
