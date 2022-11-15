subroutine calc_omega_drot_secant(vphiy,alphw,psiw,bvydw,hyydw,omega,flag_st)
  use phys_constant, only  : long
  use coordinate_grav_r, only : rg
  use def_matter_parameter, only : ome, A2DR, DRAT_A2DR, index_DRq, index_DR, &
  &                                     B2DR, DRAT_B2DR, index_DRp
  implicit none
  real(long) :: vphiy, alphw, psiw, bvydw, hyydw, omega
! -- Jome: functional of Omega (specific angular momentum in Newtonian)
  real(long) :: Jome, omegatmp, Jometmp, pmfac
  real(long) :: x_st(0:2), f_st(0:1)
  real(long) :: delta, dfdx
  real(long) :: error_st, diff_st, epsilon = 1.0d-9, small = 1.0d-30
  integer    :: ii, ic, nic, flag_st
  real(long) :: facp(5) = (/ 1.0d-1, 2.0d-1, 4.0d-1, 7.0d-1, 1.0d-0 /)

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
  nic = 40
!
  do ic = 0, nic
!
!rotlaw_type0    if(omega.gt.ome) omega = ome
!rotlaw_type2    if(omega.lt.ome) omega = ome
!
    call calc_J_utuphi(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome)
!rotlaw_type0    Jometmp = A2DR*omega*((ome/(omega+small))**index_DR - 1.0d0)
!rotlaw_type0    diff_st = Jome - Jometmp
!rotlaw_type0    error_st= dabs(diff_st/(Jome+small))
!rotlaw_type1    pmfac = -1.0d0
!rotlaw_type2    pmfac =  1.0d0
!!rotlaw_type12    Jometmp = A2DR*omega*(pmfac* &
!!rotlaw_type12    &       (omega/ome - 1.0d0)+small)**index_DR
!rotlaw_type12    Jometmp = A2DR*omega*(pmfac* &
!rotlaw_type12    &       (omega/ome - 1.0d0))
!rotlaw_type12    diff_st = Jome**(1.0/index_DR) - Jometmp
!!rotlaw_type12    error_st= dabs(diff_st/(Jometmp+small))
!rotlaw_type12    error_st= dabs((x_st(1)-x_st(0))/(x_st(1)+small))
!rotlaw_typeOJ    omegatmp = ome*(1.0d0+(Jome/(B2DR*ome))**index_DRp) &
!rotlaw_typeOJ    &             /(1.0d0+(Jome/(A2DR*ome))**(index_DRq+index_DRp))
!rotlaw_typeOJ    diff_st = omega - omegatmp
!rotlaw_typeOJ    error_st= dabs(diff_st/omega)
!rotlaw_typeJC    omegatmp = ome*(1.0d0+(Jome/(B2DR*ome))**index_DRp) &
!rotlaw_typeJC    &             *(1.0d0-(Jome/(A2DR*ome)))
!rotlaw_typeJC    diff_st = omega - omegatmp
!rotlaw_typeJC    error_st= dabs(diff_st/omega)
!!    write(6,'(1a6,1p,3e12.4)')'omeome', omega, omegatmp
!!!    write(6,'(1a6,1p,3e12.4)')'omeome', Jome, Jometmp
!!!    write(6,'(1a6,1p,3e12.4)')'omeome', diff_st, error_st
!!!    write(6,*) A2DR,omega,pmfac
!!!    write(6,*) ome,small,index_DR
!
    if (ic.ne.0) then
      x_st(0) = x_st(1)
      f_st(0) = f_st(1)
    end if
!
    if (ic.eq.0.or.error_st.gt.epsilon) then 
      x_st(1) = omega
      f_st(1) = diff_st
    else
      exit
    end if
!
    if (ic.eq.0) then
      delta = omega*0.1d0
      if (delta.eq.0.0d0) delta = 0.01d0
!rotlaw_type01      delta = - omega*0.1d0
!rotlaw_type01      if (delta.eq.0.0d0) delta = - 0.01d0
    else
      dfdx = (f_st(1) - f_st(0))/(x_st(1) - x_st(0))
      delta = - f_st(1)/dfdx
    end if
!
    ii = min0(5,ic+1)
    x_st(2) = x_st(1) + delta*facp(ii)
    omega = x_st(2)
!
!!!    write(6,'(1a6,i10)')'omseca', ic
!!!    write(6,'(1a6,1p,3e12.4)')'omseca', x_st(0),x_st(1),x_st(2)
!!!    write(6,'(1a6,1p,3e12.4)')'omseca', f_st(0),f_st(1)
!!!    write(6,'(1a6,1p,3e12.4)')'omseca', delta, dfdx, omega
!!!    write(6,*)'omseca', x_st(0),x_st(1)
!!!    write(6,*)'omseca', f_st(0),f_st(1)
!!!    write(6,*)'omseca', dfdx
  end do 
!
  flag_st = 0
  if (ic.ge.nic-2) flag_st = 1
!rotlaw_type01  if (omega.lt.0.0.or.omega.gt.ome) flag_st = 1
!rotlaw_type2  if (omega.lt.ome) flag_st = 1
!
!!write(6,*) ome,A2DR,B2DR
!!write(6,*) index_DRq,index_DRp
!!write(6,*) ic, nic, error_st
!
!!  if (ic.eq.nic) then 
!!    write(6,*) ic, omega
!!    write(6,*) 'calc_omega_drot_secant.f90_typeOJ not converged'
!!!    stop 'calc_omega_drot_secant.f90_typeOJ'
!!  end if
!
end subroutine calc_omega_drot_secant
