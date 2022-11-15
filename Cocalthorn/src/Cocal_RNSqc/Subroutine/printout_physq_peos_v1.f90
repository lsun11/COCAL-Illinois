subroutine printout_physq_peos_v1(iseq,flag_restmass)
  use phys_constant, only  :   long, pi
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, ntfpolp, sw_mass_iter
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome, radi, pinx
  use def_quantities
  implicit none
  real(long) :: fixeddlm, fixedvir, fixedvir_asymp
  real(long) :: omega
  integer, intent(in) :: iseq, flag_restmass
!
  if (iseq.eq.1) then 
    open(100,file='rnsphyseq_v1.dat',status='unknown')
  else 
    open(100,file='rnsphyseq_v1.dat',status='old', position="append")
  end if
!
  if (flag_restmass.le.1) write(100,*) '#### Solution did not converge #### '
  write(100,*) '== Sequence number == ', iseq
!
  write(100,*) '## Coordinate and Proper Radii in G = c = Msol = 1 unit ##'
  write(100,'(a22,1p,2e23.15)') ' NS radius along x  = ', &
  &                              coord_radius_x, proper_radius_x
  write(100,'(a22,1p,2e23.15)') ' NS radius along y  = ', &
  &                              coord_radius_y, proper_radius_y
  write(100,'(a22,1p,2e23.15)') ' NS radius along z  = ', &
  &                              coord_radius_z, proper_radius_z
  write(100,'(a22,1p,2e23.15)') ' Axis ratio y/x     = ', &
  &           coord_radius_y/coord_radius_x, proper_radius_y/proper_radius_x
  write(100,'(a22,1p,2e23.15)') ' Axis ratio z/x     = ', &
  &           coord_radius_z/coord_radius_x, proper_radius_z/proper_radius_x
  write(100,*) '## Coordinate and Proper Radii in [km] ##'
  write(100,'(a22,1p,2e23.15)') ' NS radius along x  = ', &
  &                              coord_radius_x_km, proper_radius_x_km
  write(100,'(a22,1p,2e23.15)') ' NS radius along y  = ', &
  &                              coord_radius_y_km, proper_radius_y_km
  write(100,'(a22,1p,2e23.15)') ' NS radius along z  = ', &
  &                              coord_radius_z_km, proper_radius_z_km
!
  write(100,*) '## Radii in Req = 1 unit ##'
  write(100,'(a22,1p,2e23.15)') ' NS radius along x  = ', rs(ntfeq,npfxzp)
  write(100,'(a22,1p,2e23.15)') ' NS radius along y  = ', rs(ntfeq,npfyzp)
  write(100,'(a22,1p,2e23.15)') ' NS radius along z  = ', rs(ntfpolp,npfxzp)
  write(100,*) '## Eccentricity in coordinate and Proper Radii ##'
  write(100,'(a22,1p,2e23.15)') ' sqrt(1-(Rz/Rx)^2)  = ', &
  &           sqrt(1.0d0-(coord_radius_z/coord_radius_x)**2), &
  &           sqrt(1.0d0-(proper_radius_z/proper_radius_x)**2)
!
  write(100,*) '## M = spherical gravitational mass (G = c = Msol = 1 unit) ##'
  write(100,'(a22,1p,2e23.15)') ' Omega M and Omega  = ', & 
  &                               ome/radi*gravmass_sph, ome/radi
  write(100,'(a22,1p,2e23.15)') ' M_ADM              = ', admmass,admmass_asymp
  write(100,'(a22,1p,2e23.15)') ' M_K Komar mass     = ', komarmass, &
  &                                                       komarmass_asymp
  write(100,'(a22,1p,2e23.15)') ' Rest mass M_0 1star= ', restmass
  write(100,'(a22,1p,2e23.15)') ' Proper mass M_p    = ', propermass
  write(100,'(a22,1p,2e23.15)') ' Angular momentum J = ', angmom, angmom_asymp
  write(100,'(a22,1p,2e23.15)') ' Spherical rest mass= ', restmass_sph
  write(100,'(a22,1p,2e23.15)') ' Spherical grav mass= ', gravmass_sph
  write(100,'(a22,1p,2e23.15)') ' Spherical M/R      = ', MoverR_sph
  write(100,'(a22,1p,2e23.15)') ' Schwarzschildradius= ', schwarz_radi_sph
  write(100,'(a22,1p,2e23.15)') ' Schwarz radi in[km]= ', schwarz_radi_sph_km
!
  write(100,'(a22,1p,2e23.15)') ' E/M = (M_ADM-M)/M  = ', &
  &                                            admmass/gravmass_sph - 1.0d0
  write(100,'(a22,1p,2e23.15)') ' J/M^2 and J/M_ADM^2= ', &
  &                     angmom/gravmass_sph**2.0d0, angmom/admmass**2.0d0
!
  write(100,'(a22,1p,2e23.15)') ' T kinetic energy   = ', T_kinene, &
  &                                                       T_kinene_omeJ 
  write(100,'(a22,1p,2e23.15)') ' P internal  energy = ', P_intene
  write(100,'(a22,1p,2e23.15)') ' W grav energy      = ', W_gravene, &
  &                                                       W_gravene_omeJ 
  write(100,'(a22,1p,2e23.15)') ' T/|W|              = ', ToverW, &
  &                                                       ToverW_omeJ
  write(100,'(a22,1p,2e23.15)') ' P/|W|              = ', PoverW
  write(100,'(a22,1p,2e23.15)') ' Virial relation    = ', Virial
  write(100,'(a22,1p,2e23.15)') ' Moment of inertia  = ', I_inertia
!
  write(100,*) '## Virial and Rest mass accuracy ##'
  fixedvir  = (admmass  -    komarmass)/admmass
  fixedvir_asymp  = (admmass_asymp  -    komarmass_asymp)/admmass_asymp
  fixeddlm  = (restmass - restmass_sph)/restmass_sph
  write(100,'(a22,1p,2e23.15)') ' 1 - M_K/M_ADM      = ', fixedvir, &
  &                                                       fixedvir_asymp
  write(100,'(a22,1p,2e23.15)') ' M_0/M_0sph - 1     = ', fixeddlm
!
  write(100,*) '## The density, pressure and p/rho at the NS center ##'
  write(100,*) '## in G = c = M = 1 unit and cgs unit          ##'
  write(100,'(a22,1p,2e23.15)') ' rho_c              = ', rho_c, rho_c_cgs
  write(100,'(a22,1p,2e23.15)') ' pre_c              = ', pre_c, pre_c_cgs
  write(100,'(a22,1p,2e23.15)') ' epsilon_c          = ', epsi_c, epsi_c_cgs
  write(100,'(a22,1p,2e23.15)') ' (p/rho)_c          = ', q_c, q_c_cgs
!
  write(100,*) '## Red and blue shift ##'
  write(100,'(a22,1p,2e23.15)') ' surface on x axis  = ', &
  &                               zrb_xp_plus, zrb_xp_minus 
  write(100,'(a22,1p,2e23.15)') ' surface on y axis  = ', &
  &                               zrb_yp_plus, zrb_yp_minus
  write(100,'(a22,1p,2e23.15)') ' surface on z axis  = ', &
  &                               zrb_zp_plus, zrb_zp_minus
!
  write(100,*) '## Ratio of enthalpy gradient at x/z and y/z ##'
  write(100,'(a22,1p,2e23.15)') ' dhdx/dhdz,dhdy/dhdz= ', &
  &                               dhdr_x/dhdr_z, dhdr_y/dhdr_z
!
  write(100,*) '## Angular velocity and frequency in cgs unit ##'
  write(100,'(a22,1p,2e23.15)') ' Omega[rad/s], f[Hz]= ', ome_cgs, ome_cgs/2.0d0/pi
!
  write(100,*) '## Quasilocal spin estimate on the surface and at nrf+1 ##'
  write(100,'(a22,1p,2e23.15)') ' Quasi local spin S = ', qua_loc_spin_surf, qua_loc_spin
  write(100,'(a22,1p,2e23.15)') ' S/M^2 and S/M_ADM^2= ',  &
  &                     qua_loc_spin/gravmass_sph**2.0d0, qua_loc_spin/admmass**2.0d0
!
  write(100,*) '## Circulation line and surface integrals at 3 planes   ##'
  write(100,'(a22,1p,2e23.15)') ' Shift xy-circul.   = ', circ_shift_xy
  write(100,'(a22,1p,2e23.15)') ' Circul at xy plane = ', circ_line_xy, circ_surf_xy
  write(100,'(a22,1p,2e23.15)') ' Circul at yz plane = ', circ_line_yz, circ_surf_yz
  write(100,'(a22,1p,2e23.15)') ' Circul at zx plane = ', circ_line_zx, circ_surf_zx
!
  close(100)
end subroutine printout_physq_peos_v1
