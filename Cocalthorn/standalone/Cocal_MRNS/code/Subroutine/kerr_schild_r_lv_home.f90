subroutine kerr_schild_r_lv_home(x,y,z,r,lvu,lvd,H)
  use phys_constant,  only : long
  use def_bh_parameter, only : m_kerr, a_kerr
  implicit none
  real(long) :: lvu(0:3), lvd(0:3), r, H
  real(long) :: x, y, z, r2xyz, a, m, r2, a2, z2
!
  m = m_kerr
  a = a_kerr
  a2 = a_kerr**2
  r2xyz = x**2 + y**2 + z**2
  z2 = z**2
!
  r  = sqrt(0.5*(r2xyz-a2 + sqrt((r2xyz-a2)**2 + 4.0d0*a2*z2)))
  r2 = r**2
!
  lvu(0) = -1.0d0
  lvu(1) = (r*x + a*y)/(r2 + a2)
  lvu(2) = (r*y - a*x)/(r2 + a2)
  lvu(3) = z/r
  lvd(0) = 1.0d0
  lvd(1) = lvu(1)
  lvd(2) = lvu(2)
  lvd(3) = lvu(3)
!
  H = 2.0d0*m*r**3/(r2**2+a2*z2) 
!
end subroutine kerr_schild_r_lv_home
