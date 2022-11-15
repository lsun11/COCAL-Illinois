module def_physq_to_array
  use phys_constant, only : long
  use def_quantities
  use def_quantities_derived
  implicit none
  real(long)  :: physq_array(100,0:100)
contains
!
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
subroutine physq_to_array(iseq,flag_restmass,num_physq)
  implicit none
  integer :: iseq, flag_restmass, num_physq, jj
!
  jj = iseq
  physq_array( 1,jj) = dble(iseq)
  physq_array( 2,jj) = dble(flag_restmass)
  physq_array( 3,jj) = coord_radius_x
  physq_array( 4,jj) = proper_radius_x
  physq_array( 5,jj) = coord_radius_y
  physq_array( 6,jj) = proper_radius_y
  physq_array( 7,jj) = coord_radius_z
  physq_array( 8,jj) = proper_radius_z
  physq_array( 9,jj) = coord_axis_ratio_yx
  physq_array(10,jj) = proper_axis_ratio_yx
  physq_array(11,jj) = coord_axis_ratio_zx
  physq_array(12,jj) = proper_axis_ratio_zx
  physq_array(13,jj) = coord_radius_normalized_x
  physq_array(14,jj) = coord_radius_normalized_y
  physq_array(15,jj) = coord_radius_normalized_z
  physq_array(16,jj) = coord_eccentricity_zx
  physq_array(17,jj) = proper_eccentricity_zx
!
  physq_array(18,jj) = omega_M
  physq_array(19,jj) = omega
  physq_array(20,jj) = admmass
  physq_array(21,jj) = komarmass
  physq_array(22,jj) = restmass
  physq_array(23,jj) = propermass
  physq_array(24,jj) = angmom
  physq_array(25,jj) = restmass_sph
  physq_array(26,jj) = gravmass_sph
  physq_array(27,jj) = MoverR_sph
  physq_array(28,jj) = schwarz_radi_sph
!
  physq_array(29,jj) = E_bind_over_M
  physq_array(30,jj) = J_over_M2
  physq_array(31,jj) = J_over_Madm2
  physq_array(32,jj) = T_kinene
  physq_array(33,jj) = W_gravene
  physq_array(34,jj) = ToverW
  physq_array(35,jj) = I_inertia
  physq_array(36,jj) = Madm_MK_virial
  physq_array(37,jj) = M0_minus_M0sph_normalized
  physq_array(38,jj) = rho_c_normalized
  physq_array(39,jj) = pre_c_normalized
  physq_array(40,jj) = epsilon_c_normalized
  physq_array(41,jj) = p_over_rho_c_normalized
  physq_array(42,jj) = zrb_xp_plus
  physq_array(43,jj) = zrb_xp_minus
  physq_array(44,jj) = zrb_yp_plus
  physq_array(45,jj) = zrb_yp_minus
  physq_array(46,jj) = zrb_zp_plus
  physq_array(47,jj) = zrb_zp_minus
  physq_array(48,jj) = omega2_M2
!
  num_physq = 48
!
end subroutine physq_to_array
!
subroutine array_to_physq(iseq,flag_restmass,num_physq)
  implicit none
  integer :: iseq, flag_restmass, num_physq, jj
!
  jj = iseq
  iseq = idnint(physq_array( 1,jj))
  flag_restmass = idnint(physq_array( 2,jj))
  coord_radius_x = physq_array( 3,jj)
  proper_radius_x = physq_array( 4,jj)
  coord_radius_y = physq_array( 5,jj)
  proper_radius_y = physq_array( 6,jj)
  coord_radius_z = physq_array( 7,jj)
  proper_radius_z = physq_array( 8,jj)
  coord_axis_ratio_yx = physq_array( 9,jj)
  proper_axis_ratio_yx = physq_array(10,jj)
  coord_axis_ratio_zx = physq_array(11,jj)
  proper_axis_ratio_zx = physq_array(12,jj)
  coord_radius_normalized_x = physq_array(13,jj)
  coord_radius_normalized_y = physq_array(14,jj)
  coord_radius_normalized_z = physq_array(15,jj)
  coord_eccentricity_zx = physq_array(16,jj)
  proper_eccentricity_zx = physq_array(17,jj)
!
  omega_M = physq_array(18,jj)
  omega = physq_array(19,jj)
  admmass = physq_array(20,jj)
  komarmass = physq_array(21,jj)
  restmass = physq_array(22,jj)
  propermass = physq_array(23,jj)
  angmom = physq_array(24,jj)
  restmass_sph = physq_array(25,jj)
  gravmass_sph = physq_array(26,jj)
  MoverR_sph = physq_array(27,jj)
  schwarz_radi_sph = physq_array(28,jj)
!
  E_bind_over_M = physq_array(29,jj)
  J_over_M2 = physq_array(30,jj)
  J_over_Madm2 = physq_array(31,jj)
  T_kinene = physq_array(32,jj)
  W_gravene = physq_array(33,jj)
  ToverW = physq_array(34,jj)
  I_inertia = physq_array(35,jj)
  Madm_MK_virial = physq_array(36,jj)
  M0_minus_M0sph_normalized = physq_array(37,jj)
  rho_c_normalized = physq_array(38,jj)
  pre_c_normalized = physq_array(39,jj)
  epsilon_c_normalized = physq_array(40,jj)
  p_over_rho_c_normalized = physq_array(41,jj)
  zrb_xp_plus = physq_array(42,jj)
  zrb_xp_minus = physq_array(43,jj)
  zrb_yp_plus = physq_array(44,jj)
  zrb_yp_minus = physq_array(45,jj)
  zrb_zp_plus = physq_array(46,jj)
  zrb_zp_minus = physq_array(47,jj)
  omega2_M2 = physq_array(48,jj)
!
  num_physq = 48
!
end subroutine array_to_physq
end module def_physq_to_array
