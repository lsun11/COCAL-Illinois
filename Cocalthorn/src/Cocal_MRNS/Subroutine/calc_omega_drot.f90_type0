subroutine calc_omega_drot(vphiy,alphw,psiw,bvydw,hyydw,omega,Jome,Jome_int)
  use phys_constant, only  : long
  use coordinate_grav_r, only : rg
  use def_matter_parameter
!  use def_vector_phi, only : vec_phi
  implicit none
  real(long), intent(in):: vphiy, alphw, psiw, bvydw, hyydw
  real(long), intent(inout):: omega
  real(long), intent(out):: Jome, Jome_int
! -- Jome: functional of Omega (specific angular momentum in Newtonian)
! -- Jome_int: integration of Jome w.r.t. omega
!
  real(long) :: sgg, gg
  real(long) :: p4a2, ovyu, ov2, ovphi, term1, term2, term3
  real(long) :: dJomedo, dterm1do, dterm2do, dterm3do
  real(long) :: omegaold, ddomega, logomega, ddlogomega 
  real(long) :: facfac, error, small, smallome = 1.0d-30
  integer    :: numero
!-------------------------------------------------------------
! Assyme axisymmetry.  
! Compute only on xz plane then copy to the other phi=const planes.
!
  small =  1.0d-30
  if (index_DR.gt.6.0d0) small = ome*0.001d0
  sgg = 0.0d0
  gg  = 0.0d0
!
  p4a2 = psiw**4/alphw**2
  numero = 1
!
  if (vphiy.lt.0.01d0*rg(1)) then
    omega = ome
  else
    logomega = dlog(omega)
!testtest
  if (index_DR.le.3.0d0) then 
  do 
    Jome = A2DR*omega*((ome/(omega+small))**index_DR - 1.0d0)
    ovyu = bvydw + omega*vphiy
! This is in fact only the \omega^y component.
!      ov2 = ovyu**2+ovxu**2+ovzu**2
    ov2   = ovyu**2*(1.0d0 + hyydw)
    ovphi = ovyu*vphiy*(1.0d0 + hyydw)
    term1 = 1.0d0 - p4a2*ov2
    term2 = p4a2*ovphi
    term3 = A2DR*omega+term2/term1
!
    dJomedo  = - A2DR + A2DR*(ome/(omega+small))**index_DR*(-index_DR+1.0d0)
    dterm1do = - p4a2*2.0d0*ovyu*vphiy
    dterm2do =   p4a2*vphiy*vphiy
    dterm3do =   A2DR + dterm2do/term1 - term2/term1**2*dterm1do
!
!testtest    if (index_DR.le.6.0d0) then 
      sgg = -(Jome*term1  - term2)
      gg  = dJomedo*term1 + Jome*dterm1do - dterm2do
      ddomega = sgg/gg
      facfac = dmin1(dble(numero)/5.0d0,1.0d0)
      omega = omega + ddomega*facfac
!testtest    else 
!test      call calc_omega_drot_bisection(vphiy,alphw,psiw,bvydw,hyydw,omega)
!!      sgg = -(dlog(Jome) + dlog(term1)  - dlog(term2))
!!      gg  = dJomedo/Jome + dterm1do/term1 - dterm2do/term2
!!      sgg = -(term1  - term2/Jome)
!!      gg  = dterm1do - dterm2do/Jome + term2*dJomedo/Jome**2
!org      sgg = -(A2DR*omega**index_DR + term2/term1*omega**(index_DR-1.0d0) & 
!org      &      -A2DR*ome**index_DR)
!org      gg  = A2DR*index_DR*omega**(index_DR-1.0d0) & 
!org      &   + dterm2do/term1*omega**(index_DR-1.0d0) & 
!org      &   - term2/term1**2*dterm1do*omega**(index_DR-1.0d0) & 
!org      &   + term2/term1*(index_DR-1.0d0)*omega**(index_DR-2.0d0)
!org      ddomega = sgg/gg
!org      facfac = dmin1(dble(numero)/5.0d0,1.0d0)
!org      omega = omega + ddomega*facfac
! take log(omega) as a variable  d(omega)/d(log(omega)) = omega
!!      sgg = -((index_DR-1.0d0)*dlog(omega) + dlog(term3) &
!!      &     - dlog(A2DR*ome**index_DR))
!!      gg  = (index_DR-1.0d0)/omega + (dterm3do/term3)
!log      sgg = -((index_DR-1.0d0)*dlog(omega) + dlog(term3) &
!log      &     - dlog(A2DR*ome**index_DR))
!log      gg  = (index_DR-1.0d0) + (dterm3do/term3)*omega
!log      ddlogomega = sgg/gg
!log      facfac = dmin1(dble(numero)/5.0d0,1.0d0)
!log      logomega = logomega + ddlogomega*facfac
!log      ddomega  = dexp(logomega) - omega
!log      omega    = dexp(logomega)
!testtest    end if
!
    if (omega.gt.ome) omega = ome
    if (omega.le.0.0d0) omega = smallome
!!    error = dabs(ddomega/ome)
    error = dabs(ddomega/omega)
    numero = numero + 1
    if (numero>1000) then
write(6,*) vphiy, ome
write(6,*) ome, p4a2, bvydw
write(6,*) Jome, ovyu, omega
      write(6,*)' numero = ', numero, '   error =',error
    end if
!      if (iter<10.and.numero>10) exit
    if (numero>1010) stop 'rotation law'
    if (dabs(error)<1.0d-14) exit
!write(6,*) Jome, ovyu, omega
!!write(6,*) Jome1, dJome1do
!!write(6,*) Jome2, dJome2do
!write(6,*) omega,ddomega,sgg
!!pause
  end do
    else
      call calc_omega_drot_bisection(vphiy,alphw,psiw,bvydw,hyydw,omega)
    end if
  end if
!!  end if
!
  Jome = A2DR*omega*((ome/omega)**index_DR - 1.00d0)
  if (index_DR.ne.2.0d0) then
    Jome_int = A2DR*omega**2*((ome/omega)**index_DR/(2.0d0-index_DR) - 0.5d0) &
    &        - index_DR*A2DR*ome**2/(4.0d0-2.0d0*index_DR)
  else ! v-constant law
    Jome_int = - A2DR*ome**2*dlog(ome/omega) + 0.5d0*A2DR*(ome**2-omega**2)
  end if
!
end subroutine calc_omega_drot
