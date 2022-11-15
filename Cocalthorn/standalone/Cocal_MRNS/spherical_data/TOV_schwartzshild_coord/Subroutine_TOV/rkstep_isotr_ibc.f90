subroutine rkstep_isotr(neq, func, x, h, y0, yn, ak, w, rfac)
  implicit none
  real(8), parameter   :: a2 = 0.5d0, a3 = a2, &
                       &  b2 = 0.5d0, b3 = b2, &
                       &  c1 = 1.0d0 / 6, c2 = 1.0d0 / 3, c3 = c2, c4 = c1
  integer, intent(in)    :: neq
  real(8), intent(in)    :: rfac
  real(8), intent(inout) :: x, h, y0(neq), ak(neq)
  real(8), intent(out)   :: yn(neq), w(neq)
  integer :: i
!
  call func(x, y0, ak, rfac)
  yn(1:neq) = y0(1:neq) + h * c1 * ak(1:neq)
  w(1:neq) = y0(1:neq) + h * b2 * ak(1:neq)
  call func(x + a2 * h, w, ak, rfac)
  do i = 1, neq
    yn(i) = yn(i) + h * c2 * ak(i)
  end do
  w(1:neq) = y0(1:neq) + h * b3 * ak(1:neq)
  call func(x + a3 * h, w, ak, rfac)
  do i = 1, neq
    yn(i) = yn(i) + h * c3 * ak(i)
  end do
  w(1:neq) = y0(1:neq) + h * ak(1:neq)
  call func(x + h, w, ak, rfac)
  do i = 1, neq
    yn(i) = yn(i) + h * c4 * ak(i)
  end do
  return
end subroutine rkstep_isotr
