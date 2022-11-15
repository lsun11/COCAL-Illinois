subroutine peos_h2qprho(h,q,pre,rho,ened)
  use def_peos_parameter   !rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(inout)  :: h
  real(8), intent(out)    :: q, pre, rho, ened
  real(8), external       :: lagint_4th, lagint_2nd
  real(8) :: x4(4), f4(4), x2(2), f2(2)
  integer              :: i, i0, ii
  integer, save        :: iphase
  real(8)              :: qq, hin, qin, abin, abct
  real(8)              :: fac1, fac2, fack, small = 1.0d-16
  character(2) :: char1
!
!  if (dabs(h-1.0d0) < 1.0d-14) then
!    q=small;  rho=small;  pre=q*rho;  ened=rho*h-pre
!    q=0.0d0;  rho=0.0d0;  pre=0.0d0;  ened=0.0d0
!    return
!  end if

  if(h<=hi(0))  then
    iphase=1

    hin  = hi(0) 
    qin  = qi(0)
    abin = gamma_crust
    abct = kappa_crust
    fac1 = 1.0d0/(abin - 1.0d0)
    fac2 = abin/(abin - 1.0d0)
    fack = abct**(-fac1)
    q  = (h - hin)/fac2 + qin
    qq = q
    if (h <= 1.0d0 .or. qq==0.0d0) qq = small/fac2
    pre = fack*qq**fac2
    rho = fack*qq**fac1
    ened = rho*h - pre

!    write(6,*) "hin,qin,abin,abct: ", hin,qin, abin,abct
!    write(6,*) "q,qq,rho,pre: ", q,qq,rho,pre
    return
  else
    ii=iphase
    if(h.lt.hi(iphase)) then
      do i = ii, 1, -1
        if( (h.ge.hi(i-1)) .and. (h.le.hi(i))) then
          iphase = i
          char1 = 'lt'
            exit
        end if
      end do
    else
      do i = ii, nphase
        if( (h.ge.hi(i-1)) .and. (h.le.hi(i))) then
          iphase = i
          char1 = 'ge'
          exit
        end if
      end do
    end if
  end if

!  do i = 1, nphase
!    if( (h.ge.hi(i-1)) .and. (h.le.hi(i))) then
!      iphase = i
!      exit
!    end if
!  end do
!   write(6,*) "h2qprho: ii, iphase=", ii, iphase, char1
!
  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = hi(i0:i0+3)
  f4(1:4) = qi(i0:i0+3)
  q       = lagint_4th(x4,f4,h)

!write(6,'(a10,i5,1p,2e23.15)') "i0,h = ", i0, h
!write(6,'(a5,1p,5e23.15)') "x4:", (x4(i),i=1,4)
!write(6,'(a5,1p,5e23.15)') "q4:", (f4(i),i=1,4)
!write(6,*) " q=", q

  f4(1:4) = prei(i0:i0+3)
  pre     = lagint_4th(x4,f4,h)

!write(6,*) "pre4:", (f4(i),i=1,4)
!write(6,*) " pre=", pre

  f4(1:4) = rhoi(i0:i0+3)
  rho     = lagint_4th(x4,f4,h)

!write(6,'(a5,1p,5e23.15)') "rho4:", (f4(i),i=1,4)
!write(6,'(a5,1p,5e23.15)') " rho=", rho

  f4(1:4) = enei(i0:i0+3)
  ened    = lagint_4th(x4,f4,h)

!write(6,*) "ene4:", (f4(i),i=1,4)
!write(6,*) " ene=", ened

!
end subroutine peos_h2qprho
