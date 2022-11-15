! computes the spherical bessel functions j_l(x), y_l(x),
! with y_l = n_l the spherical Neumann function 
! The routine uses the Numerical Recipes routine bessjy 
! input: x, l
! output: sj = j_l(x), sy=y_l(x)
! -----------------------------------------------------------
subroutine sphbess(x,l,sj,sy)
  implicit none
  real(8) :: x, sj, sy
  integer :: l
  real(8) :: rtpi2 = 1.253314154753822d0
  real(8) :: xnu, xrj, xry, xrjp, xryp
!
  xnu = dble(l) + 0.5d0
  call bessjy(x,xnu,xrj,xry,xrjp,xryp)
  sj = xrj*rtpi2/sqrt(x)
  sy = xry*rtpi2/sqrt(x)
end subroutine sphbess
