subroutine rk(neq, func, x0, xe, n, y0, yn, work)
!*********************************************************************
!    subroutine rk numerically integrates a system of neq
!    first order ordinary differential equations of the form
!            dy(i)/dx = f(x, y(1), ..., y(neq)),
!    by the classical runge-kutta formula.
!
!    parameters
!  === input ===
!    (1) neq: number of equations to be integrated
!    (2) func: subroutine func(x, y, f) to evaluate derivatives
!               f(i)=dy(i)/dx
!    (3) x0: initial value of independent variable
!    (4) xe: output point at which the solution is desired
!    (5) n: number of divisions
!       the interval (x0, xe) is divided into n subintervals
!       with the length (xe-x0)/n and in each subinterval
!       the classical runge-kutta formula is used.
!    (6) y0(i) (i=1, .., neq): initial value at x0
!  === output ===
!    (7) yn(i) (i=1, .., neq): approximate solution at xe
!  === other ===
!    (8) work(): two-dimentional array (size=(neq, 2)) to be
!                used inside rk
!    copyright: m. sugihara, november 15, 1989, v. 1
!*********************************************************************
  implicit none
  external func
  integer, intent(in)    :: neq, n
  real(8), intent(inout) :: x0, y0(neq), work(neq, 2)
  real(8), intent(in)    :: xe
  real(8), intent(out)   :: yn(neq)
  real(8) :: h
  integer :: i
  h = (xe - x0) / n
  do i = 1, n
    call rkstep(neq, func, x0, h, y0, yn, work(1, 1), work(1, 2))
    x0 = x0 + h
    y0(1:neq) = yn(1:neq)
  end do
  x0 = xe
  return
end subroutine rk
