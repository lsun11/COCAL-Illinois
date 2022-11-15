subroutine IO_input_physq_plot(iseq,flag_restmass)
  use def_quantities
  use def_quantities_derived
  implicit none
  integer, intent(out) :: iseq, flag_restmass
!
  if (iseq.eq.1) then 
    open(210,file='rnsphyplot_all_input.dat',status='old')
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
! 9  10 ' Axis ratio y/x     = ', coord_axis_ratio_yx, proper_axis_ratio_yx
! 11 12 ' Axis ratio z/x     = ', coord_axis_ratio_zx, proper_axis_ratio_zx
!    '## Radii in Req = 1 unit ##'
! 13 ' NS radius along x  = ', coord_radius_normalized_x
! 14 ' NS radius along y  = ', coord_radius_normalized_y
! 15 ' NS radius along z  = ', coord_radius_normalized_z
!    '## Eccentricity in coordinate and Proper Radii ##'
! 16 17 ' sqrt(1-(Rz/Rx)^2)  = ', coord_eccentricity_zx,proper_eccentricity_zx
!
!    '## M = spherical gravitational mass (K = 1 unit) ##'
! 18 19 ' Omega M and Omega  = ', omega_M, omega
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
! 29    ' E/M = (M_ADM-M)/M  = ', E_bind_over_M
! 30 31 ' J/M^2 and J/M_ADM^2= ', J_over_M2, J_over_Madm2
!
! 32 ' T kinetic energy   = ', T_kinene
! 33 ' W grav energy      = ', W_gravene
! 34 ' T/|W|              = ', ToverW
! 35 ' Moment of inertia  = ', I_inertia
!
! '## Virial and Rest mass accuracy ##'
! 36  ' 1 - M_K/M_ADM      = ', Madm_MK_virial
! 37  ' M_0/M_0sph - 1     = ', M0_minus_M0sph_normalized
!
! '## Maximums of the density, pressure and p/rho ##'
! 38 ' rho_max            = ', rho_c_normalized
! 39 ' pre_max            = ', pre_c_normalized
! 40 ' epsilon_max        = ', epsilon_c_normalized
! 41 ' (p/rho)_max        = ', p_over_rho_c_normalized
!
! '## Red and blue shift ##'
! 42  43 ' surface on x axis  = ', zrb_xp_plus, zrb_xp_minus 
! 44, 45 ' surface on y axis  = ', zrb_yp_plus, zrb_yp_minus
! 46, 47 ' surface on z axis  = ', zrb_zp_plus, zrb_zp_minus
!
! 48 ' (Omega M)**2 '   omega2_M2
!
!
read(210,'(2i5,1p,48e23.15)',ERR=9999,END=9999) &
iseq, flag_restmass, &
coord_radius_x, proper_radius_x, &
coord_radius_y, proper_radius_y, &
coord_radius_z, proper_radius_z, &
coord_axis_ratio_yx, proper_axis_ratio_yx, &
coord_axis_ratio_zx, proper_axis_ratio_zx, &
coord_radius_normalized_x, &
coord_radius_normalized_y, &
coord_radius_normalized_z, &
coord_eccentricity_zx, proper_eccentricity_zx, &
!
omega_M, omega, &
admmass, komarmass, restmass, propermass, angmom, restmass_sph, & 
gravmass_sph, MoverR_sph, schwarz_radi_sph, & 
!
E_bind_over_M, J_over_M2, J_over_Madm2, &
T_kinene, W_gravene, ToverW, I_inertia, &
Madm_MK_virial, M0_minus_M0sph_normalized, &
rho_c_normalized, pre_c_normalized, &
epsilon_c_normalized, p_over_rho_c_normalized, &
zrb_xp_plus,zrb_xp_minus,zrb_yp_plus,zrb_yp_minus,zrb_zp_plus,zrb_zp_minus, &
omega2_M2
!
return
9999 continue
flag_restmass = 9999
return
!
end subroutine IO_input_physq_plot
