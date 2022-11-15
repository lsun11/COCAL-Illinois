subroutine printout_TOV_quantities_isotr_tidal
  use def_matter_parameter, only : emdc, pinx
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  implicit none

  real(8) :: yri, xs,xscgs,  xi,xicgs, y,co2,co3,co5,deno
!
  yr = yn_iti(9)/yn_iti(1)*(1.0d0-(1.0d0-2.0d0*yn_iti(1)/yn_iti(9))**0.5)  !Correct
  yri = 1.0d0 + yn_iti(1)/2.0d0/xe   !Wrong since xe not correct radius in isotropic coordinates
!
  xs = yn_iti(9)  ! Correct radius in Schwarzschild coordinates
  xi = 0.5d0*(xs*(1.0d0+(1.0d0-2.0d0*yn_iti(1)/xs)**0.5) - yn_iti(1)) !Correct radius in isotropic coord
!
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a55,1es23.15)') "1. Radius in isotropic coord xe                      = ", xe
  write(6,'(a55,1es23.15)') "2. Radius in Schwarzschild coord                yn(9)= ", yn_iti(9)
  write(6,'(a55,1es23.15)') "3. Mass                                         yn(1)= ", yn_iti(1)
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a55,1es23.15)') "4. ri = 0.5d0*(r*(1.0d0+(1.0d0-2.0d0*m/r)**0.5) - m) = ", xi
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a34,1es23.15)') "5. Psi from r Schwarzschild  2. = ", yr
  write(6,'(a34,1es23.15)') "6. Psi from r isotropic (xe) 1. = ", yri
  write(6,'(a34,1es23.15)') "7. Psi calculated         yn(5) = ", yn_iti(5)
  write(6,'(a34,1es23.15)') "8. Psi from isotropic ri 4.     = ",  1.0d0 + yn_iti(1)/2.0d0/xi
  write(6,*)'------------------------------------------------------------------------------'
  write(6,'(a36,1es23.15)') "9. Alpha from r Schwarzschild  2. = ", (1.0d0 - 2.0d0*yn_iti(1)/yn_iti(9))**0.5
  write(6,'(a36,1es23.15)') "10.Alpha from r isotropic  1.     = ",          &
    &                        (1.0d0 - yn_iti(1)/2.0d0/xe)/(1.0d0 + yn_iti(1)/2.0d0/xe)
  write(6,'(a36,1es23.15)') "11.Alpha calculated         yn(8) = ", yn_iti(8)
  write(6,*)'------------------------------------------------------------------------------'
!
  adm =  yn_iti(5)*yn_iti(6)/yr 
  compa = yn_iti(1)/yn_iti(9)   ! compactness in Schwarzschild coord
!
  y  = xs*yn_iti(12)/yn_iti(11)

  co2=compa*compa
  co3=co2*compa
  co5=co3*co2

  deno = 2.0d0*compa*(6.0d0-3.0d0*y+3.0d0*compa*(5.0d0*y-8.0d0))  &
   &   + 4.0d0*co3*(13.0d0-11.0d0*y+compa*(3.0d0*y-2.0d0)+2.0d0*co2*(1.0d0+y)) &
   &   + 3.0d0*(1.0d0-2.0d0*compa)**2*(2.0d0-y+2.0d0*compa*(y-1.0d0))*dlog(1.0d0-2.0d0*compa)

  k2 = 8.0d0*co5/5.0d0*(1.0d0-2.0d0*compa)**2*(2.0d0+2.0d0*compa*(y-1.0d0)-y)/deno

!
  call peos_h2qprho_tidal(hini, qini, pre, rho0, ene, dpde)
!
  emdc = qini
  restmass_sph = yn_iti(3)
  gravmass_sph = yn_iti(1)
  propermass_sph = yn_iti(4)
  MoverR_sph  = compa
  schwarz_radi_sph = yn_iti(9)
!
  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  radicgs = radi*g*solmas*1.0d-5/c**2
  xscgs = xs*g*solmas*1.0d-5/c**2
  xicgs = xi*g*solmas*1.0d-5/c**2
!
  open(3,file='ovphy_plot.dat',status='unknown')
  write(3 , '(20es14.6)') compa, emdc, rho0, rhocgs, pre, ene, &
  &       xi, xicgs, restmass_sph, propermass_sph, adm, gravmass_sph, yn_iti(6), &
  &       xs, xscgs, y
  close(3)
!
  open(20,file='ovphy.dat',status='unknown')
  write(20,'(a50,1es23.15)') ' Spherical M/R                                  = ', compa
  write(20,'(a50,1es23.15)') ' p/rho                                          = ', emdc
!
  write(20,'(a50,1es23.15)') ' rho (G=c=Msol=1 unit)                          = ', rho0
  write(20,'(a50,1es23.15)') ' rho (cgs)                                      = ', rhocgs
!
  write(20,'(a50,1es23.15)') ' pre(G=c=Msol=1 unit)                           = ', pre
  write(20,'(a50,1es23.15)') ' pre(cgs)                                       = ', precgs
!
  write(20,'(a50,1es23.15)') ' epsilon(G=c=Msol=1 unit)                       = ',ene
  write(20,'(a50,1es23.15)') ' epsilon(cgs)                                   = ',epsiloncgs
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
  write(20,'(a50,1es23.15)') ' adm yn_iti(6)                                  = ', yn_iti(6)
  write(20,'(a50,1es23.15)') ' y = R*dH/dr(R)/H(R)                            = ', y
  write(20,'(a50,1es23.15)') ' k2=3*G*lambda/(2R**5)                          = ', k2
  write(20,'(a50)') ''

  write(20,'(a50)') ''
  write(20,'(a67,a72,a50)')  "$\ M_{ADM}[M_\odot]\ $ & $\ M_0[M_\odot]\ $ & $\ M_p[M_\odot]\ $ & ",   &
    &                "$\ R[km]\ $ & $\quad\ M_{ADM}/R\quad $ & $\quad \log(\rho_0)_c\quad $ & ", &
    &                "$\quad \log P_c\quad $ &  $\quad \log\rho_c\quad$ "

  write(20,'(f6.3,a3,f6.3,a3,f6.3,a3,f7.3,a3,f6.3,a3,f7.3,a3,f7.3,a3,f7.3)') adm, &
   &   ' & ', restmass_sph, ' & ', propermass_sph, ' & ', xscgs, ' & ', compa, ' & ', &
   &  dlog10(rhocgs), ' & ', dlog10(precgs), ' & ', dlog10(epsiloncgs)

  close(20)
!
end subroutine printout_TOV_quantities_isotr_tidal
