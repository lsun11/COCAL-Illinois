function chebev(a,b,c,m,x)
  implicit none
  integer :: m, j
  real(8) :: a, b, x, c(m)
  real(8) :: d, dd, sv, y, y2
  real(8) :: chebev
  if ((x-a)*(x-b).gt.0.0d0) stop 'x not in range in chebev'
  d  = 0.0d0
  dd = 0.0d0
  y  =(2.0d0*x-a-b)/(b-a)
  y2 = 2.0d0*y
  do j = m, 2, -1
    sv = d
    d = y2*d-dd+c(j)
    dd = sv
  end do
  chebev = y*d - dd + 0.5d0*c(1)
  return
end function chebev
!
