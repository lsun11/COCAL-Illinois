subroutine calc_omega_drot(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome,Jome_int)
  use phys_constant, only  : long
  use coordinate_grav_r, only : rg
  use def_matter_parameter
!  use def_vector_phi, only : vec_phi
  implicit none
  real(long), intent(in):: vphiy, alphw, psiw, bvydw, hyydw
  real(long), intent(inout):: omega, Jome, Jome_int
! -- Jome: functional of Omega (specific angular momentum in Newtonian)
! -- Jome_int: integration of Jome w.r.t. omega
!
  real(long) :: omega_bak, pmfac
  integer    :: flag_st
!-------------------------------------------------------------
! Assume axisymmetry.  
! Compute only on xz plane then copy to the other phi=const planes.
!
  omega_bak = omega
!
  flag_st = 1
!rotlaw_OJ2nd  call calc_omega_drot_secant(vphiy,alphw,psiw,bvydw,hyydw,omega,flag_st)
!rotlaw_OJ4th  call calc_omega_drot_secant(vphiy,alphw,psiw,bvydw,hyydw,omega,flag_st)
!rotlaw_OJjco  call calc_omega_drot_secant(vphiy,alphw,psiw,bvydw,hyydw,omega,flag_st)
!
!  call calc_omega_drot_secant(vphiy,alphw,psiw,bvydw,hyydw,omega,flag_st)
  if (flag_st.eq.1) then
    omega = omega_bak
    call calc_omega_drot_bisection(vphiy,alphw,psiw,bvydw,hyydw,omega)
  end if
!
!rotlaw_type0  Jome = A2DR*omega*((ome/omega)**index_DR - 1.00d0)
!rotlaw_type1  pmfac = -1.0d0
!rotlaw_type1  Jome = A2DR*omega*(pmfac*(omega/ome - 1.0d0))**index_DR
!rotlaw_type2  pmfac =  1.0d0
!rotlaw_type2  Jome = A2DR*omega*(pmfac*(omega/ome - 1.0d0))**index_DR
!rotlaw_OJ2nd  call calc_J_utuphi(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome)
!rotlaw_OJ4th  call calc_J_utuphi(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome)
!rotlaw_OJjco  call calc_J_utuphi(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome)
!
!rotlaw_type0  call calc_Jome_int('type0',omega,Jome,Jome_int)
!rotlaw_type1  call calc_Jome_int('type1',omega,Jome,Jome_int)
!rotlaw_type2  call calc_Jome_int('type2',omega,Jome,Jome_int)
!rotlaw_OJ2nd  call calc_Jome_int('OJ2nd',omega,Jome,Jome_int)
!rotlaw_OJ4th  call calc_Jome_int('OJ4th',omega,Jome,Jome_int)
!rotlaw_OJjco  call calc_Jome_int('OJjco',omega,Jome,Jome_int)
!
end subroutine calc_omega_drot
