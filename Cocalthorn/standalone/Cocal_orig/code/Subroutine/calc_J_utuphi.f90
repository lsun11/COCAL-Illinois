subroutine calc_J_utuphi(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome)
  use phys_constant, only  : long
  implicit none
  real(long), intent(in) :: vphiy, alphw, psiw, bvydw, hyydw, omega
  real(long), intent(out):: Jome
  real(long) :: p4a2, ovyu, ov2, ovphi, term1, term2
!
  p4a2 = psiw**4/alphw**2
  ovyu = bvydw + omega*vphiy
! This is in fact only the \omega^y component.
!      ov2 = ovyu**2+ovxu**2+ovzu**2
  ov2   = ovyu**2*(1.0d0 + hyydw)
  ovphi = ovyu*vphiy*(1.0d0 + hyydw)
  term1 = 1.0d0 - p4a2*ov2
  term2 = p4a2*ovphi
  Jome = term2/term1
!
end subroutine calc_J_utuphi
