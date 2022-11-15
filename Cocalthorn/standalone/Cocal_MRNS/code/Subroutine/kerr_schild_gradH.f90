subroutine kerr_schild_gradH(x,y,z,r,H,gradH)
  use phys_constant,  only : long
  use def_bh_parameter, only : m_kerr, a_kerr
  implicit none
  real(long) :: x, y, z, r, H, gradH(3)
  real(long) :: m, m2, a, a2, H2, H3, H3m2, H2mr32
!
  m = m_kerr; m2 = m_kerr**2
  a = a_kerr; a2 = a_kerr**2
  H2= H**2  ; H3 = H**3
  H3m2   = H3/m2
  H2mr32 = 3.0d0*H2/(2.0d0*m*r)
!
  gradH(1) = -x*H3m2 + x*H2mr32
  gradH(2) = -y*H3m2 + y*H2mr32
  gradH(3) = -z*H3m2 + z*H2mr32 &
  &        +(0.5d0*m - H*r)*a2*H2*z/(m2*r**3)
!
end subroutine kerr_schild_gradH
