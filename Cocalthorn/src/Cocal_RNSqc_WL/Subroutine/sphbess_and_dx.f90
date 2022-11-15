! computes the spherical bessel functions j_l(x), y_l(x),
! with y_l = n_l the spherical Neumann function 
! The routine uses the Numerical Recipes routine bessjy 
! input: x, l
! output: sj = j_l(x), sy=y_l(x)
! -----------------------------------------------------------
subroutine sphbess_and_dx(x,l,sj,sy,dxsj,dxsy)
  implicit none
  real(8) :: x, sj, sy, dxsj, dxsy
  integer :: l
  real(8) :: rtpi2 = 1.253314154753822d0
  real(8) :: xnu, xrj, xry, xrjp, xryp
  real(8) :: oneover2x, rtpi2overx
!
  xnu = dble(l) + 0.5d0
  call bessjy(x,xnu,xrj,xry,xrjp,xryp)
  oneover2x = 0.5d0/x
  rtpi2overx = rtpi2/sqrt(x)
  sj = xrj*rtpi2overx
  sy = xry*rtpi2overx
  dxsj = rtpi2overx*(xrjp - oneover2x*xrj)
  dxsy = rtpi2overx*(xryp - oneover2x*xry)
end subroutine sphbess_and_dx
