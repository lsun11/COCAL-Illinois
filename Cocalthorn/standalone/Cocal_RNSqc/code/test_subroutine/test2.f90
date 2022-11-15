!
  implicit none
!
  real(8) :: pi = 3.14159265358979d+0
  real(8) :: x(5), y(5)
  integer :: ii, jj
!
  do ii = 1, 5
    x(ii) = -1.0 + 0.5*dble(ii-1)
    y(ii) = -1.0 + 0.5*dble(ii-1)
  end do
  do jj = 1, 5
  do ii = 1, 5
  write(6,'(5es11.3)') x(ii),y(jj),atan2(y(jj),x(ii))/pi, &
 &                     dmod(2.0*pi + atan2(y(jj),x(ii)),2.0*pi)/pi
  end do
  end do
!
  end
