subroutine peos_lookup(qp,qpar,iphase)
  use phys_constant             !nnteos
  use def_peos_parameter        !rhoi,qi,hi,nphase
  implicit none
!
  real(8), intent(in)  :: qp, qpar(0:nnteos)
  real(8)              :: det
  integer, intent(out) :: iphase
  integer              :: ii, i
!
! --  Monotonically increasing qpar is assumed.
!
!  ii=iphase
!  if(qp.lt.qpar(iphase)) then
!    do i = ii, 1, -1
!      if( (qp.ge.qpar(i-1)) .and. (qp.le.qpar(i)))  then
!        iphase = i
!        exit
!      end if
!    end do
!  else
!    do i = ii, nphase
!      if( (qp.ge.qpar(i-1)) .and. (qp.le.qpar(i)))  then
!        iphase = i
!        exit
!      end if
!    end do
!  endif

  iphase = 1
  do ii = 1, nphase
    det = (qp-qpar(ii))*(qp-qpar(ii-1))
    if (det <= 0.0d0) then
      iphase = ii
      exit
    end if
  end do

end subroutine peos_lookup
