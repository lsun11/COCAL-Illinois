subroutine calc_omega_drot_WL(vphix,vphiy,al,ps,bx,by,bz,  &
      &       hxx,hxy,hxz,hyy,hyz,hzz,omega,Jome,Jome_int)
  use phys_constant, only  : long
  use coordinate_grav_r, only : rg
  use def_matter_parameter
  implicit none
  real(long), intent(in):: vphix, vphiy, al, ps, bx, by, bz, &
    &                      hxx, hxy, hxz, hyy, hyz, hzz
  real(long), intent(inout):: omega
  real(long), intent(out):: Jome, Jome_int
! -- Jome: functional of Omega (specific angular momentum in Newtonian)
! -- Jome_int: integration of Jome w.r.t. omega
!
  real(long) :: sgg, gg, gxx,gxy,gxz,gyy,gyz,gzz
  real(long) :: p4a2, ovxu,ovyu,ovzu, ov2, ovphi, phi2, term1, term2
  real(long) :: dJomedo, dterm1do, dterm2do
  real(long) :: omegaold, ddomega, logomega, ddlogomega 
  real(long) :: facfac, error, small, smallome = 1.0d-30
  integer    :: numero
!-------------------------------------------------------------
! Assume axisymmetry.  
! Compute only on xz plane then copy to the other phi=const planes.
!
  small =  1.0d-30
  if (index_DR.gt.6.0d0) small = ome*0.001d0
  sgg = 0.0d0
  gg  = 0.0d0
!
  p4a2 = ps**4/al**2
  numero = 1
!
  if (index_DR.le.3.0d0) then 
    do 
!
!     **************** The rotation law ***********************************************
      Jome = A2DR*omega*((ome/(omega+small))**index_DR - 1.0d0)
      dJomedo  = - A2DR + A2DR*(ome/(omega+small))**index_DR*(-index_DR+1.0d0)
!     *********************************************************************************
!
      gxx=1.0d0+hxx; gxy=hxy; gxz=hxz
      gyy=1.0d0+hyy; gyz=hyz
      gzz=1.0d0+hzz
!
      ovxu = bx + omega*vphix
      ovyu = by + omega*vphiy
      ovzu = bz 
!
      ov2   = gxx*ovxu*ovxu + gxy*2.0d0*ovxu*ovyu + gxz*2.0d0*ovxu*ovzu &
       &    + gyy*ovyu*ovyu + gyz*2.0d0*ovyu*ovzu + gzz*ovzu*ovzu
      ovphi = gxx*ovxu*vphix + gxy*(ovxu*vphiy+ovyu*vphix) + gxz*ovzu*vphix &
       &    + gyy*ovyu*vphiy + gyz*ovzu*vphiy 
      phi2  = gxx*vphix*vphix + gxy*2.0d0*vphix*vphiy + gyy*vphiy*vphiy

      term1 = 1.0d0 - p4a2*ov2
      term2 = p4a2*ovphi
!
      dterm1do = - p4a2*2.0d0*ovphi
      dterm2do =   p4a2*phi2
!
      sgg = -(Jome*term1  - term2)
      gg  = dJomedo*term1 + Jome*dterm1do - dterm2do
      ddomega = sgg/gg
      facfac = dmin1(dble(numero)/5.0d0,1.0d0)
      omega = omega + ddomega*facfac

!aaa    if (omega.gt.ome) omega = ome
!aaa      if (omega.le.0.0d0) omega = smallome
      if (omega.le.0.0d0) omega = -omega
      error = dabs(ddomega/omega)
      numero = numero + 1

!        write(6,*) "--------------------------------------------"
!        write(6,*) ps, al
!        write(6,*) bx,by,bz
!        write(6,*) hxx,hxy,hxz
!        write(6,*) hyy,hyz,hzz
!        write(6,*) ov2, ovphi, phi2
!        write(6,*) ddomega, ome, omega, Jome
!        write(6,*)' numero = ', numero, '   error =',error
 
      if (numero>50) then
        write(6,*) ps, al
        write(6,*) bx,by,bz
        write(6,*) hxx,hxy,hxz
        write(6,*) hyy,hyz,hzz
        write(6,*) ov2, ovphi, phi2
        write(6,*) ddomega, ome, omega, Jome
        write(6,*)' numero = ', numero, '   error =',error
      end if
!      if (iter<10.and.numero>10) exit
      if (numero>52) stop 'rotation law'
      if (dabs(error)<1.0d-14) exit
    end do
  else
    call calc_omega_drot_bisection(vphiy,al,ps,by,hyy,omega)
  end if
!
  Jome = A2DR*omega*((ome/omega)**index_DR - 1.00d0)
  if (index_DR.ne.2.0d0) then
    Jome_int = A2DR*omega**2*((ome/omega)**index_DR/(2.0d0-index_DR) - 0.5d0) &
    &        - index_DR*A2DR*ome**2/(4.0d0-2.0d0*index_DR)
  else ! v-constant law
    Jome_int = - A2DR*ome**2*dlog(ome/omega) + 0.5d0*A2DR*(ome**2-omega**2)
  end if
!
end subroutine calc_omega_drot_WL
