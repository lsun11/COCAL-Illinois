subroutine printout_TOV_quantities_isotr_mit
  use def_matter_parameter, only : emdc, pinx
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  use def_qeos_parameter
  implicit none

  real(8) :: yri, xs,xscgs,  xi,xicgs
!
  yr = yn_iso(9)/yn_iso(1)*(1.0d0-(1.0d0-2.0d0*yn_iso(1)/yn_iso(9))**0.5)  !Correct
  yri = 1.0d0 + yn_iso(1)/2.0d0/xe   !Wrong since xe not correct radius in isotropic coordinates
!
  xs = yn_iso(9)  ! Correct radius in Schwarzschild coordinates
  xi = 0.5d0*(xs*(1.0d0+(1.0d0-2.0d0*yn_iso(1)/xs)**0.5) - yn_iso(1)) !Correct radius in isotropic coord
!
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a55,1es23.15)') "1. Radius in isotropic coord xe                      = ", xe
  write(6,'(a55,1es23.15)') "2. Radius in Schwarzschild coord                yn(9)= ", yn_iso(9)
  write(6,'(a55,1es23.15)') "3. Mass                                         yn(1)= ", yn_iso(1)
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a55,1es23.15)') "4. ri = 0.5d0*(r*(1.0d0+(1.0d0-2.0d0*m/r)**0.5) - m) = ", xi
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a34,1es23.15)') "5. Psi from r Schwarzschild  2. = ", yr
  write(6,'(a34,1es23.15)') "6. Psi from r isotropic (xe) 1. = ", yri
  write(6,'(a34,1es23.15)') "7. Psi calculated         yn(5) = ", yn_iso(5)
  write(6,'(a34,1es23.15)') "8. Psi from isotropic ri 4.     = ",  1.0d0 + yn_iso(1)/2.0d0/xi
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a36,1es23.15)') "9. Alpha from r Schwarzschild  2. = ", (1.0d0 - 2.0d0*yn_iso(1)/yn_iso(9))**0.5
  write(6,'(a36,1es23.15)') "10.Alpha from r isotropic  1.     = ",          &
    &                        (1.0d0 - yn_iso(1)/2.0d0/xe)/(1.0d0 + yn_iso(1)/2.0d0/xe)
  write(6,'(a36,1es23.15)') "11.Alpha calculated         yn(8) = ", yn_iso(8)
  write(6,*)'------------------------------------------------------------------------------'
!
  adm =  yn_iso(5)*yn_iso(6)/yr 
  compa = yn_iso(1)/yn_iso(9)   ! compactness in Schwarzschild coord
!
  ene = eneini_gcm1
  call mit_ene2prho(ene, pre, rho0)
!  qini = pre/ene
!
!  emdc = qini
  restmass_sph = yn_iso(3)
  gravmass_sph = yn_iso(1)
  propermass_sph = yn_iso(4)
  MoverR_sph  = compa
  schwarz_radi_sph = yn_iso(9)
!
  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  radicgs = radi*g*solmas*1.0d-5/c**2
  xscgs = xs*g*solmas*1.0d-5/c**2
  xicgs = xi*g*solmas*1.0d-5/c**2
  rho0 = eneini_cgs/(c**2)
!
  open(3,file='ovphy_plot.dat',status='unknown')
  write(3 , '(14es14.6)') compa, eneini_cgs, ene,enesurf_cgs,rho0, pre,  &
  &       xi, xicgs, restmass_sph, propermass_sph, adm, gravmass_sph, yn_iso(6), xs
  close(3)
!
  open(20,file='ovphy.dat',status='unknown')
  write(20,'(a50,1es23.15)') ' Spherical M/R                                  = ', compa
  write(20,'(a50,1es23.15)') ' central energy density (cgs)                   = ', eneini_cgs
!
  write(20,'(a50,1es23.15)') ' central energy density (G=c=Msol=1 unit)       = ', ene
  write(20,'(a50,1es23.15)') ' surface energy density (cgs)                   = ', enesurf_cgs
  write(20,'(a50,1es23.15)') ' surface energy density (G=c=Msol=1 unit)       = ', enesurf_gcm1
!
  write(20,'(a50,1es23.15)') ' pre(G=c=Msol=1 unit)                           = ', pre
  write(20,'(a50,1es23.15)') ' pre(cgs)                                       = ', precgs
!
!
  write(20,'(a50,1es23.15)') ' Calculated isotropic Schradius(G=c=Msol=1 unit)= ',radi
  write(20,'(a50,1es23.15)') ' Correct isotropic radius(G=c=Msol=1 unit)      = ',xi
  write(20,'(a50,1es23.15)') ' Correct isotropic Schwarzschildradius(km)      = ',xicgs
!
  write(20,'(a50,1es23.15)') ' Schwarzschildradius(G=c=Msol=1 unit) Schwarzsc = ',xs
  write(20,'(a50,1es23.15)') ' Schwarzschildradius(km) Schwarzschild coord    = ',xscgs
!
  write(20,'(a50,1es23.15)') ' Spherical rest mass                            = ',restmass_sph
  write(20,'(a50,1es23.15)') ' Proper mass M_p                                = ',propermass_sph
  write(20,'(a50,1es23.15)') ' M_ADM                                          = ', adm 
  write(20,'(a50,1es23.15)') ' Gravitational mass M                           = ', gravmass_sph
  write(20,'(a50,1es23.15)') ' adm yn_iso(6)                                  = ', yn_iso(6)
  write(20,'(a50)') ''
  close(20)
!
end subroutine printout_TOV_quantities_isotr_mit
