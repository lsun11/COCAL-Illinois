!
  implicit none
!
  real(8) :: x(4), f(4), y(1:10,1:20)
  real(8) :: rrff
  integer :: ii, ir0, jj
!
  do ii = 1, 10
    x(ii) = dble(ii)
  end do
  do jj = 1, 10
  do ii = 1, 10
    y(ii,jj) = 10. + dble(ii)
  end do
  end do
!
  ir0 = 3
  f(1:4) = 2.0*y(ir0:ir0+3,9)
  write(6,'(1p,5e11.3)') f(1:4)
  end
