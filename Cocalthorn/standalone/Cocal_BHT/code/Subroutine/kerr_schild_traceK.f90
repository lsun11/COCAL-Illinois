subroutine kerr_schild_traceK(x,y,z,traceK)
  use phys_constant,  only : long
  use def_bh_parameter,  only : m_kerr
  implicit none
  real(long) :: traceK
  real(long) :: x, y, z, r, lvu(0:3), lvd(0:3), H, gradH(3)
  real(long) :: H2, Hp1, H2H1, H21pH, lvugradH, m, a
!
  call kerr_schild_r_lv_home(x,y,z,r,lvu,lvd,H)
  call kerr_schild_gradH(x,y,z,r,H,gradH)
  m = m_kerr
  H2  = H**2 ; Hp1 = 1.0d0+H 
  H2H1  = 0.5d0*(2.0d0+H)/Hp1**1.5d0
  H21pH = H**2/sqrt(Hp1)
  lvugradH = lvu(1)*gradH(1)+lvu(2)*gradH(2)+lvu(3)*gradH(3)
!
  traceK = H2H1*lvugradH + H21pH/m
!
end subroutine kerr_schild_traceK
