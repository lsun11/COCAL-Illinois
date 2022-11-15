subroutine kerr_schild_transverse_part(x,y,z,dgamma)
  use phys_constant,  only : long
  use def_bh_parameter,  only : m_kerr
  implicit none
  real(long) :: dgamma(3)
  real(long) :: x, y, z, r, lvu(0:3), lvd(0:3), H, gradH(3)
  real(long) :: Hp1, Hp1m23, H3H1, san, lvugradH, m, H2m
!
  san = 1.0d0/3.0d0 
  call kerr_schild_r_lv_home(x,y,z,r,lvu,lvd,H)
  call kerr_schild_gradH(x,y,z,r,H,gradH)
  m = m_kerr ; H2m = H**2/m
  Hp1      = 1.0d0 + H
  Hp1m23   = Hp1**(-2.0d0/3.0d0)
  H3H1     =(3.0d0 + H)/Hp1
  lvugradH = lvu(1)*gradH(1)+lvu(2)*gradH(2)+lvu(3)*gradH(3)
!
  dgamma(1:3) = san*Hp1m23*(gradH(1:3) - H3H1*lvugradH*lvu(1:3)) &
  &           - H2m*Hp1m23*lvu(1:3)
!
end subroutine kerr_schild_transverse_part
