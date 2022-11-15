subroutine printout_physq_plot(iseq,flag_restmass)
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, ntfpolp, sw_mass_iter
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome, radi, pinx
  use def_quantities
  implicit none
  real(8) :: fixeddlm, fixedvir
  integer, intent(in) :: iseq, flag_restmass
!
  if (iseq.eq.1) then 
    open(200,file='rnsphyplot_all.dat',status='unknown')
!    open(210,file='rnsphyplot.dat',status='unknown')
  else
    open(200,file='rnsphyplot_all.dat',status='old', position="append")
  end if
! 
! Reference for the data 
! Slot # 
! 1  '== Sequence number == ', iseq
! 2  '#### Solution did not converge #### ', flag_restmass
!    '## Coordinate and Proper Radii in K = 1 unit ##'
! 3, 4  ' NS radius along x  = ',  coord_radius_x, proper_radius_x
! 5, 6  ' NS radius along y  = ',  coord_radius_y, proper_radius_y
! 7, 8  ' NS radius along z  = ',  coord_radius_z, proper_radius_z
! 9  ' Axis ratio y/x     = ', coord_radius_y/coord_radius_x
! 10                           proper_radius_y/proper_radius_x
! 11 ' Axis ratio z/x     = ', coord_radius_z/coord_radius_x
! 12                           proper_radius_z/proper_radius_x
!    '## Radii in Req = 1 unit ##'
! 13 ' NS radius along x  = ', rs(ntfeq,npfxzp)
! 14 ' NS radius along y  = ', rs(ntfeq,npfyzp)
! 15 ' NS radius along z  = ', rs(ntfpolp,npfxzp)
!    '## Eccentricity in coordinate and Proper Radii ##'
! 16 ' sqrt(1-(Rz/Rx)^2)  = ', sqrt(1.0d0-(coord_radius_z/coord_radius_x)**2)
! 17                           sqrt(1.0d0-(proper_radius_z/proper_radius_x)**2)
!
!    '## M = spherical gravitational mass (K = 1 unit) ##'
! 18 19 ' Omega M and Omega  = ', ome/radi*gravmass_sph, ome/radi
! 20    ' M_ADM              = ', admmass
! 21    ' M_K Komar mass     = ', komarmass
! 22    ' Rest mass M_0 1star= ', restmass
! 23    ' Proper mass M_p    = ', propermass
! 24    ' Angular momentum J = ', angmom
! 25    ' Spherical rest mass= ', restmass_sph
! 26    ' Spherical grav mass= ', gravmass_sph
! 27    ' Spherical M/R      = ', MoverR_sph
! 28    ' Schwarzschildradius= ', schwarz_radi_sph
!
! 29 ' E/M = (M_ADM-M)/M  = ', admmass/gravmass_sph - 1.0d0
! 30 ' J/M^2 and J/M_ADM^2= ', angmom/gravmass_sph**2.0d0
! 31                           angmom/admmass**2.0d0
!
! 32 ' T kinetic energy   = ', T_kinene
! 33 ' W grav energy      = ', W_gravene
! 34 ' T/|W|              = ', ToverW
! 35 ' Moment of inertia  = ', I_inertia
!
! '## Virial and Rest mass accuracy ##'
! 36  ' 1 - M_K/M_ADM      = ', fixedvir  = (admmass  -    komarmass)/admmass
! 37  ' M_0/M_0sph - 1     = ', fixeddlm  = (restmass - restmass_sph)
!                                                       /restmass_sphfixeddlm
! '## Maximums of the density, pressure and p/rho ##'
! 38 ' rho_max            = ', rho_max
! 39 ' pre_max            = ', pre_max
! 40 ' epsilon_max        = ', epsi_max
! 41 ' (p/rho)_max        = ', q_max
!
! '## Red and blue shift ##'
! 42  43 ' surface on x axis  = ', zrb_xp_plus, zrb_xp_minus 
! 44, 45 ' surface on y axis  = ', zrb_yp_plus, zrb_yp_minus
! 46, 47 ' surface on z axis  = ', zrb_zp_plus, zrb_zp_minus
!
! 48 ' (Omega M)**2 '   (ome/radi*gravmass_sph)**2
!
! 49 ' Omega[rad/s] = ', ome_cgs
! 50 ' Quasi local spin S = ', qua_loc_spin
! 51, 52 ' S/M^2 and S/M_ADM^2= ',  qua_loc_spin/gravmass_sph**2.0d0, qua_loc_spin/admmass**2.0d0

  fixedvir  = (admmass  -    komarmass)/admmass
  fixeddlm  = (restmass - restmass_sph)/restmass_sph
!
write(200,'(2i5,1p,52e23.15)') &
iseq, flag_restmass, &
coord_radius_x, proper_radius_x, &
coord_radius_y, proper_radius_y, &
coord_radius_z, proper_radius_z, &
coord_radius_y/coord_radius_x, proper_radius_y/proper_radius_x,  &
coord_radius_z/coord_radius_x, proper_radius_z/proper_radius_x,  &
rs(ntfeq,npfxzp), rs(ntfeq,npfyzp), rs(ntfpolp,npfxzp), &
sqrt(1.0d0-(coord_radius_z/coord_radius_x)**2), &
sqrt(1.0d0-(proper_radius_z/proper_radius_x)**2), &
!
ome/radi*gravmass_sph, ome/radi, & 
admmass, komarmass, restmass, propermass, angmom, restmass_sph, & 
gravmass_sph, MoverR_sph, schwarz_radi_sph, & 
!
admmass/gravmass_sph - 1.0d0, &
angmom/gravmass_sph**2.0d0, angmom/admmass**2.0d0, &                     
T_kinene, W_gravene, ToverW, I_inertia, fixedvir, fixeddlm, &
rho_max, pre_max, epsi_max, q_max, &
zrb_xp_plus,zrb_xp_minus,zrb_yp_plus,zrb_yp_minus,zrb_zp_plus,zrb_zp_minus, &
(ome/radi*gravmass_sph)**2, &
ome_cgs, &
qua_loc_spin, qua_loc_spin/gravmass_sph**2.0d0, qua_loc_spin/admmass**2.0d0

!
  close(200)
!
end subroutine printout_physq_plot
