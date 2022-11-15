subroutine calc_omega_drot_Newton(vphiy,alphw,psiw,bvydw,hyydw,omega,flag_st)
  use phys_constant, only  : long
  use coordinate_grav_r, only : rg
  use def_matter_parameter, only : ome, A2DR, DRAT_A2DR, index_DR
  implicit none
  real(long) :: vphiy, alphw, psiw, bvydw, hyydw, omega
! -- Jome: functional of Omega (specific angular momentum in Newtonian)
  real(long) :: Jome, pmfac
  real(long) :: sgg, gg
  real(long) :: p4a2, ovyu, ov2, ovphi, term1, term2, term3
  real(long) :: dJomedo, dterm1do, dterm2do
  real(long) :: ddomega
  real(long) :: facfac, error, small = 1.0d-30
  integer    :: ic, nic, flag_st
  real(long) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)
!
!-------------------------------------------------------------
! Assume axisymmetry.  
!
  if (ome  .eq.0.0d0) ome   = 0.01d0 ! For 1D initial
  if (omega.eq.0.0d0) omega = 0.01d0 ! For 1D initial
!
  if (vphiy.lt.0.001d0*rg(1)) then
    omega = ome
    return
  end if
!
!rotlaw_type1  pmfac = -1.0d0
!rotlaw_type2  pmfac =  1.0d0
!
  nic = 30
!
  do ic = 0, nic
!
    Jome = A2DR*omega*((ome/(omega+small))**index_DR - 1.0d0)
    ovyu = bvydw + omega*vphiy
! This is in fact only the \omega^y component.
!      ov2 = ovyu**2+ovxu**2+ovzu**2
    ov2   = ovyu**2*(1.0d0 + hyydw)
    ovphi = ovyu*vphiy*(1.0d0 + hyydw)
    term1 = 1.0d0 - p4a2*ov2
    term2 = p4a2*ovphi
!
    dJomedo  = - A2DR + A2DR*(ome/(omega+small))**index_DR*(-index_DR+1.0d0)
    dterm1do = - p4a2*2.0d0*ovyu*vphiy
    dterm2do =   p4a2*vphiy*vphiy
!
    sgg = -(Jome*term1  - term2)
    gg  = dJomedo*term1 + Jome*dterm1do - dterm2do
    ddomega = sgg/gg
    facfac = facp(min(ic+1,5))
    omega = omega + ddomega*facfac
!
    if (omega.gt.ome) omega = ome
    if (omega.le.0.0d0) omega = small
    error = dabs(ddomega/omega)
    if (dabs(error)<1.0d-14) exit
!write(6,*) Jome, ovyu, omega
!write(6,*) dJomedo
!!write(6,*) Jome1, dJome1do
!!write(6,*) Jome2, dJome2do
!write(6,*) omega,ddomega,sgg
!!pause
  end do
!
  flag_st = 0
  if (ic.ge.nic-2) flag_st = 1
!!write(6,*) ome,A2DR,B2DR
!!write(6,*) index_DRq,index_DRp
!!write(6,*) ic, nic, error_st
!
!!  if (ic.eq.nic) then 
!!    write(6,*) ic, omega
!!    write(6,*) 'calc_omega_drot_secant.f90_typeOJ not converged'
!!!    stop 'calc_omega_drot_Newton.f90_type1type2'
!!  end if
!
end subroutine calc_omega_drot_Newton
