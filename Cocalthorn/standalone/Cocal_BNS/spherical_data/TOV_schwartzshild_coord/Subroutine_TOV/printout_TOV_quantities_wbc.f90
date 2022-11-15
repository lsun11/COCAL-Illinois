subroutine printout_TOV_quantities_wbc
  use def_matter_parameter, only : emdc, pinx
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  implicit none
!
  yr = xe/yn_wbc(1)*(1.0d0-(1.0d0-2.0d0*yn_wbc(1)/xe)**0.5)
  adm =  yn_wbc(5)*yn_wbc(6)/yr
  compa = yn_wbc(1)/xe
!
!  write(6,*) xe
!  write(6,*) yn_wbc(1)
!  write(6,*) yn_wbc(5)
!  write(6,*) yr

  call peos_h2qprho(hini, qini, pre, rho0, ene)
!
  emdc = qini
  restmass_sph = yn_wbc(3)
  gravmass_sph = yn_wbc(1)
  propermass_sph = yn_wbc(4)
  MoverR_sph  = compa
  schwarz_radi_sph = xe
!
  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  radicgs = radi*g*solmas*1.0d-5/c**2
!
  open(3,file='ovphy_plot.dat',status='unknown')
  write(3 , '(13es14.6)') compa, emdc, rho0, rhocgs, pre, ene, &
  &       radi, radicgs, restmass_sph, propermass_sph, adm, gravmass_sph, yn_wbc(6)
  close(3)
!
  open(20,file='ovphy.dat',status='unknown')
  write(20,'(a39,1es23.15)') ' Spherical M/R                       = ', compa
  write(20,'(a39,1es23.15)') ' p/rho                               = ', emdc
!
  write(20,'(a39,1es23.15)') ' rho (G=c=Msol=1 unit)               = ', rho0
  write(20,'(a39,1es23.15)') ' rho (cgs)                           = ', rhocgs
!
  write(20,'(a39,1es23.15)') ' pre(G=c=Msol=1 unit)                = ', pre
  write(20,'(a39,1es23.15)') ' pre(cgs)                            = ', precgs
!
  write(20,'(a39,1es23.15)') ' epsilon(G=c=Msol=1 unit)            = ',ene
  write(20,'(a39,1es23.15)') ' epsilon(cgs)                        = ',epsiloncgs
!
  write(20,'(a39,1es23.15)') ' Schwarzschildradius(G=c=Msol=1 unit)= ',radi
  write(20,'(a39,1es23.15)') ' Schwarzschildradius(km)             = ',radicgs
!
  write(20,'(a39,1es23.15)') ' Spherical rest mass                 = ',restmass_sph
  write(20,'(a39,1es23.15)') ' Proper mass M_p                     = ',propermass_sph
  write(20,'(a39,1es23.15)') ' M_ADM                               = ', adm 
  write(20,'(a39,1es23.15)') ' Gravitational mass M                = ', gravmass_sph
  write(20,'(a39,1es23.15)') ' adm yn_wbc(6)                       = ', yn_wbc(6)
  close(20)
!
end subroutine printout_TOV_quantities_wbc
