subroutine calc_omega_drot_bisection_WL(vphix,vphiy,al,ps,bx,by,bz,  &
      &                                  hxx,hxy,hxz,hyy,hyz,hzz,omega)
  use phys_constant, only  : long
  use coordinate_grav_r, only : rg
  use def_matter_parameter
!  use def_vector_phi, only : vec_phi
  implicit none
  real(long), intent(in):: vphix, vphiy, al, ps, bx, by, bz, &
    &                      hxx, hxy, hxz, hyy, hyz, hzz
  real(long), intent(inout):: omega
  real(long) :: Jome, omega0, omega1, omegam
  real(long) :: sgg0, sgg1, sggm, sgg0m, sgg1m
  real(long) :: gxx,gxy,gxz,gyy,gyz,gzz
  real(long) :: p4a2, ovxu,ovyu,ovzu, ov2, ovphi, phi2, term1, term2, term3
  real(long) :: ddomega
  real(long) :: facfac, error, small, smallome = 1.0d-30
  integer    :: i, numero
!-------------------------------------------------------------
! Assume axisymmetry.  
!
!  if (ome.eq.0.0d0) ome = 0.01d0 ! For 1D initial
!  small = ome*0.001
  p4a2 = ps**4/al**2
  numero = 1
!
  omega0 = smallome
  omega1 = ome
  do 
    omegam = (omega0 + omega1)*0.5d0
    do i = 0, 2
      if (i.eq.0) omega = omega0
      if (i.eq.1) omega = omega1
      if (i.eq.2) omega = omegam
!
      Jome = A2DR*omega*((ome/omega)**index_DR - 1.0d0)


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
      if (i.eq.0) sgg0 = Jome*term1  - term2
      if (i.eq.1) sgg1 = Jome*term1  - term2
      if (i.eq.2) sggm = Jome*term1  - term2
    end do
    sgg0m = sgg0*sggm
    sgg1m = sgg1*sggm
!
    ddomega = omega1 - omega0
    omega = omegam
    if (sggm.eq.0.0d0) then 
      ddomega = 0.0d0
    else if (sgg0m.gt.0.0d0.and.sgg1m.lt.0.0d0) then
      omega0 = omegam
    else if (sgg0m.lt.0.0d0.and.sgg1m.gt.0.0d0) then
      omega1 = omegam
    else
      write(6,*) omega0, omegam, omega1
      write(6,*) sgg0, sggm, sgg1
      stop 'rotation law'
    end if
!
!    if (omega.gt.ome) omega = ome
!    if (omega.le.0.0d0) omega = small
!!    error = dabs(ddomega/ome)
    error = dabs(ddomega/omega)
    numero = numero + 1
    if (numero>100) then
      write(6,*) omega0, omegam, omega1
      write(6,*) vphiy, ome, ddomega
      write(6,*) ome, p4a2
      write(6,*) Jome, ovyu, omega
      write(6,*)' numero = ', numero, '   error =',error
    end if
!      if (iter<10.and.numero>10) exit
    if (numero>105) stop 'rotation law'
    if (dabs(error)<1.0d-14) exit
!write(6,*) Jome, ovyu, omega
!!write(6,*) Jome1, dJome1do
!!write(6,*) Jome2, dJome2do
!write(6,*) omega,ddomega,sgg
!!pause
  end do
!
end subroutine calc_omega_drot_bisection_WL
