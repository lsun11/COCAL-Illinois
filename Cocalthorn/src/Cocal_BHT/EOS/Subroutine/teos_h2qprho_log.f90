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
  real(8) :: small = 1.0d-17
  real(8) :: lh,lq,lpre,lrho,lened
  character(2) :: char1
!
  if (h<1.0d0) then
    q=h*small;  rho=h*small**2.8d0;  pre=q*rho;  ened=rho*h-pre
!  if (dabs(h-1.0d0) < 1.0d-14) then
!    q=small;  rho=small;  pre=q*rho;  ened=rho*h-pre
    return
  end if

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

!  do i = 1, nphase
!    if( (h.ge.hi(i-1)) .and. (h.le.hi(i))) then
!      iphase = i
!      exit
!    end if
!  end do
!   write(6,*) "h2qprho: ii, iphase=", ii, iphase, char1
!
  lh = dlog10(h)

  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = dlog10(hi(i0:i0+3))
  f4(1:4) = dlog10(qi(i0:i0+3))
  lq      = lagint_4th(x4,f4,lh)
  q       = 10.0d0**lq

!write(6,'(a5,1p,5e23.15)') "x4:", (x4(i),i=1,4)
!write(6,'(a5,1p,5e23.15)') "q4:", (f4(i),i=1,4)
!write(6,*) " q=", q

  f4(1:4) = dlog10(prei(i0:i0+3))
  lpre    = lagint_4th(x4,f4,lh)
  pre     = 10.0d0**lpre
!write(6,*) "pre4:", (f4(i),i=1,4)
!write(6,*) " pre=", q

  f4(1:4) = dlog10(rhoi(i0:i0+3))
  lrho    = lagint_4th(x4,f4,lh)
  rho     = 10.0d0**lrho 
!write(6,*) "rho4:", (f4(i),i=1,4)
!write(6,*) " rho=", q

  f4(1:4) = dlog10(enei(i0:i0+3))
  lened   = lagint_4th(x4,f4,lh)
  ened    = 10.0d0**lened
!write(6,*) "ene4:", (f4(i),i=1,4)
!write(6,*) " ene=", q

!
end subroutine peos_h2qprho
