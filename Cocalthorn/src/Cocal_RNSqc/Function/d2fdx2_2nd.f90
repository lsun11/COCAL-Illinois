function d2fdx2_2nd(xg,fnc,rv)
  implicit none
  real(8) :: xg(1:4),fnc(1:4), rv, d2fdx2_2nd
  real(8) :: dr01, dr02, dr03, dr12, dr13, &
  &          dr23, dr10, dr20, dr21, dr30, dr31, dr32, &
  &          wer0, wer1, wer2, wer3, &
  &          rrv0, rrv1, rrv2, rrv3
  integer :: ir0, ir1, ir2, ir3
!
  ir0 = 1
  ir1 = 2
  ir2 = 3
  ir3 = 4
  dr01 = xg(ir0) - xg(ir1)
  dr02 = xg(ir0) - xg(ir2)
  dr03 = xg(ir0) - xg(ir3)
  dr12 = xg(ir1) - xg(ir2)
  dr13 = xg(ir1) - xg(ir3)
  dr23 = xg(ir2) - xg(ir3)
  dr10 = - dr01
  dr20 = - dr02
  dr21 = - dr12
  dr30 = - dr03
  dr31 = - dr13
  dr32 = - dr23
  rrv0 = rv - xg(ir0)
  rrv1 = rv - xg(ir1)
  rrv2 = rv - xg(ir2)
  rrv3 = rv - xg(ir3)
  wer0 = 2.0d0*(rrv2*rrv3 + rrv1*rrv3 + rrv1*rrv2) &
     &          /(dr01*dr02*dr03)
  wer1 = 2.0d0*(rrv2*rrv3 + rrv0*rrv3 + rrv0*rrv2) &
     &          /(dr10*dr12*dr13)
  wer2 = 2.0d0*(rrv1*rrv3 + rrv0*rrv3 + rrv0*rrv1) &
     &          /(dr20*dr21*dr23)
  wer3 = 2.0d0*(rrv1*rrv2 + rrv0*rrv2 + rrv0*rrv1) &
     &          /(dr30*dr31*dr32)
!
  d2fdx2_2nd = wer0*fnc(ir0) + wer1*fnc(ir1) &
     &       + wer2*fnc(ir2) + wer3*fnc(ir3)
!
end function d2fdx2_2nd
