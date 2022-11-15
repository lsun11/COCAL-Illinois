subroutine printout_TOV_quantities
  use def_matter_parameter, only : emdc, pinx
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  implicit none
!
  yr = xe/yn(1)*(1.0d0-(1.0d0-2.0d0*yn(1)/xe)**0.5)
  adm =  yn(5)*yn(6)/yr
  compa = yn(1)/xe
!
  call peos_h2qprho(hini, qini, pre, rho0, ene)
!
  emdc = qini
  restmass_sph = yn(3)
  gravmass_sph = yn(1)
  propermass_sph = yn(4)
  MoverR_sph  = compa
  schwarz_radi_sph = xe
!
  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  radicgs = radi*g*solmas*1.0d-5/c**2
!
  open(3,file='ovphy_plot.dat',status='unknown')
  write(3 , '(11es14.6)') compa, emdc, rho0, rhocgs, pre, ene, &
  &                       radi, radicgs, restmass_sph, propermass_sph,adm
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
  close(20)
!
      open(14, file='rnspar_add.dat', status='unknown') 
!
      write(14, '(2es14.6, a14)') emdc, pinx, &
        & '  : emdc, pinx'
      write(14, '(2es14.6, a30)') restmass_sph, gravmass_sph, &
        & '  : restmass_sph, gravmass_sph'
      write(14, '(2es14.6, a49)') MoverR_sph, schwarz_radi_sph, &
        & '  : MoverR_sph,  schwarz_radi_sph  (G=c=M=1 unit)'
      write(6, '(/, /, 2es14.6, a14)') emdc, pinx, &
        & '  : emdc, pinx'
      write(6, '(2es14.6, a30)') restmass_sph, gravmass_sph, &
        & '  : restmass_sph, gravmass_sph'
      write(6, '(2es14.6, a49)') MoverR_sph, schwarz_radi_sph, &
        & '  : MoverR_sph,  schwarz_radi_sph  (G=c=M=1 unit)'
!
      close(14)
!
end subroutine printout_TOV_quantities
