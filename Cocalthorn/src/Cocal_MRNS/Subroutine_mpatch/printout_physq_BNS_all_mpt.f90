subroutine printout_physq_BNS_all_mpt(filename)
  use phys_constant,  only : nmpt, g, c, solmas
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, ntfpolp, sw_mass_iter
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome, radi, pinx
  use def_quantities
  use def_quantities_mpt
  use def_matter_parameter_mpt
  use grid_parameter_mpt
  use def_binary_parameter
  use def_binary_parameter_mpt
  implicit none
  real(8) :: MM = solmas, LL = g*solmas/c**2, TT = g*solmas/c**3
  real(8) :: Msph_tot
  real(8) :: fixeddlm, fixedvir
  integer :: i,ia,ib
  character(30) :: filename
!
  Msph_tot = def_quantities_real_(28,1) + def_quantities_real_(28,2)

  open(100,file=filename,status='unknown')

  write(100,'(a62)') '## If units are not mentioned, G=c=Msol=1 is implied        ##'
  write(100,'(a62)') '## M is total spherical gravitational mass (G=c=Msol=1)     ##'
  write(100,'(a62)') '## R0 is the rescaling factor of all coordinate patches     ##'

  write(100,'(a72)') '========================================================================'

  write(100,'(a24,1p,1e23.15)') 'Omega*R0              = ', def_matter_param_real_(3,1)

  write(100,'(a24,1p,1e23.15)') 'Omega                 = ', &
                              & def_matter_param_real_(3,1)/def_matter_param_real_(5,1)

  write(100,'(a24,1p,1e23.15)') 'Omega M_ADM           = ', &
                              & def_matter_param_real_(3,1)/def_matter_param_real_(5,1)*def_quantities_real_(7,nmpt) 

  write(100,'(a24,1p,1e23.15)') 'Omega [rad/sec]       = ', def_quantities_real_(78,1)

  write(100,'(a24,1p,1e23.15)') 'R0                    = ', def_matter_param_real_(5,1)

  write(100,'(a24,1p,2e23.15)') 'r_surf                = ', surf_param_real_(1,1), surf_param_real_(1,2)

  write(100,'(a24,1p,2e23.15)') 'Radius Rx [km]        = ', &
                              & surf_param_real_(1,1)*def_matter_param_real_(5,1)*LL*1.0d-5, &
                              & surf_param_real_(1,2)*def_matter_param_real_(5,2)*LL*1.0d-5

  write(100,'(a24,1p,2e23.15)') 'Separation            = ', sepa_(1)   

  write(100,'(a24,1p,2e23.15)') 'Separation*R0         = ', sepa_(1)*def_matter_param_real_(5,1)

  write(100,'(a24,1p,2e23.15)') 'Separation*R0 [km]    = ', sepa_(1)*def_matter_param_real_(5,1)*LL*1.0d-5

  write(100,'(a24,1p,2e23.15)') 'Distance to CM (d)    = ', dis_(1), dis_(2)

  write(100,'(a24,1p,2e23.15)') 'd*R0                  = ', dis_(1)*def_matter_param_real_(5,1), &
                                                         &  dis_(2)*def_matter_param_real_(5,2)

  write(100,'(a24,1p,2e23.15)') 'd*R0 [km]             = ', dis_(1)*def_matter_param_real_(5,1)*LL*1.0d-5, &
                                                         &  dis_(2)*def_matter_param_real_(5,2)*LL*1.0d-5

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15,a30)') 'admmass               = ', def_quantities_real_(1,1), &
                                  & def_quantities_real_(1,2), '     nc: Not valid for COCP1,2'
  write(100,'(a24,1p,2e23.15)')     'komarmass             = ', def_quantities_real_(2,1), &
                                  & def_quantities_real_(2,2)
  write(100,'(a24,1p,2e23.15,a30)') 'komarmass_nc          = ', def_quantities_real_(3,1), &
                                  & def_quantities_real_(3,2), '     nc: Not valid for COCP1,2'
  write(100,'(a24,1p,2e23.15)')     'restmass              = ', def_quantities_real_(4,1), def_quantities_real_(4,2)
  write(100,'(a24,1p,2e23.15)')     'propermass            = ', def_quantities_real_(5,1), def_quantities_real_(5,2)
  write(100,'(a24,1p,2e23.15)')     'angmom                = ', def_quantities_real_(6,1), def_quantities_real_(6,2)
  write(100,'(a24,1p,2e23.15)')     'qua_loc_spin          = ', def_quantities_real_(79,1), def_quantities_real_(79,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'admmass_asymp         = ', def_quantities_real_(7,nmpt)
  write(100,'(a24,1p,2e23.15)') 'komarmass_asymp       = ', def_quantities_real_(8,nmpt)
  write(100,'(a24,1p,2e23.15)') 'angmom_asymp          = ', def_quantities_real_(9,nmpt)
  write(100,'(a24,1p,2e23.15)') 'admmom_asymp(1)       = ', def_quantities_real_(10,nmpt)
  write(100,'(a24,1p,2e23.15)') 'admmom_asymp(2)       = ', def_quantities_real_(11,nmpt)
  write(100,'(a24,1p,2e23.15)') 'admmom_asymp(3)       = ', def_quantities_real_(12,nmpt)
  write(100,'(a24,1p,2e23.15)') 'E/M = (M_ADM-M)/M     = ', def_quantities_real_(7,nmpt)/Msph_tot - 1.0d0
  write(100,'(a24,1p,2e23.15)') 'J/M^2                 = ', def_quantities_real_(9,nmpt)/Msph_tot**2.0d0
  write(100,'(a24,1p,2e23.15)') 'J/M_ADM^2             = ', &
                              & def_quantities_real_(9,nmpt)/def_quantities_real_(7,nmpt)**2.0d0
  write(100,'(a24,1p,2e23.15)') '1 - M_K/M_ADM         = ', &
                              & 1.0d0 - def_quantities_real_(8,nmpt)/def_quantities_real_(7,nmpt) 

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'T_kinene              = ', def_quantities_real_(13,1), def_quantities_real_(13,2)
  write(100,'(a24,1p,2e23.15)') 'W_gravene             = ', def_quantities_real_(14,1), def_quantities_real_(14,2)
  write(100,'(a24,1p,2e23.15)') 'P_intene              = ', def_quantities_real_(15,1), def_quantities_real_(15,2)
  write(100,'(a24,1p,2e23.15)') 'M_emfene              = ', def_quantities_real_(16,1), def_quantities_real_(16,2)
  write(100,'(a24,1p,2e23.15)') 'M_torBene             = ', def_quantities_real_(17,1), def_quantities_real_(17,2)
  write(100,'(a24,1p,2e23.15)') 'M_polBene             = ', def_quantities_real_(18,1), def_quantities_real_(18,2)
  write(100,'(a24,1p,2e23.15)') 'M_eleEene             = ', def_quantities_real_(19,1), def_quantities_real_(19,2)

  write(100,'(a24,1p,2e23.15)') 'Virial                = ', def_quantities_real_(20,1), def_quantities_real_(20,2)

  write(100,'(a24,1p,2e23.15)') 'ToverW                = ', def_quantities_real_(21,1), def_quantities_real_(21,2)
  write(100,'(a24,1p,2e23.15)') 'PoverW                = ', def_quantities_real_(22,1), def_quantities_real_(22,2)
  write(100,'(a24,1p,2e23.15)') 'MoverW                = ', def_quantities_real_(23,1), def_quantities_real_(23,2)
  write(100,'(a24,1p,2e23.15)') 'MtorBoverW            = ', def_quantities_real_(24,1), def_quantities_real_(24,2)
  write(100,'(a24,1p,2e23.15)') 'MpolBoverW            = ', def_quantities_real_(25,1), def_quantities_real_(25,2)
  write(100,'(a24,1p,2e23.15)') 'MeleEoverW            = ', def_quantities_real_(26,1), def_quantities_real_(26,2)

  write(100,'(a24,1p,2e23.15)') 'I_inertia             = ', def_quantities_real_(27,1), def_quantities_real_(27,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'gravmass_sph          = ', def_quantities_real_(28,1), def_quantities_real_(28,2)
  write(100,'(a24,1p,2e23.15)') 'restmass_sph          = ', def_quantities_real_(29,1), def_quantities_real_(29,2)
  write(100,'(a24,1p,2e23.15)') 'propermass_sph        = ', def_quantities_real_(30,1), def_quantities_real_(30,2)
  write(100,'(a24,1p,2e23.15)') 'MoverR_sph            = ', def_quantities_real_(31,1), def_quantities_real_(31,2)
  write(100,'(a24,1p,2e23.15)') 'M_0/M_0sph - 1        = ', def_quantities_real_(4,1)/def_quantities_real_(29,1) - 1.0d0, &
                                                          & def_quantities_real_(4,2)/def_quantities_real_(29,2) - 1.0d0

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'schwarz_radi_sph      = ', def_quantities_real_(32,1), def_quantities_real_(32,2)
  write(100,'(a24,1p,2e23.15)') 'coord_radius_x        = ', def_quantities_real_(33,1), def_quantities_real_(33,2)
  write(100,'(a24,1p,2e23.15)') 'coord_radius_y        = ', def_quantities_real_(34,1), def_quantities_real_(34,2)
  write(100,'(a24,1p,2e23.15)') 'coord_radius_z        = ', def_quantities_real_(35,1), def_quantities_real_(35,2)
  write(100,'(a24,1p,2e23.15)') 'Coord axis ratio y/x  = ', def_quantities_real_(34,1)/def_quantities_real_(33,1),  &
                                                         &  def_quantities_real_(34,2)/def_quantities_real_(33,2)      
  write(100,'(a24,1p,2e23.15)') 'Coord axis ratio z/x  = ', def_quantities_real_(35,1)/def_quantities_real_(33,1),  &
                                                         &  def_quantities_real_(35,2)/def_quantities_real_(33,2)
  write(100,'(a24,1p,2e23.15)') 'sqrt(1-(Ry/Rx)^2)     = ', &
                              & sqrt(1.0d0-(def_quantities_real_(34,1)/def_quantities_real_(33,1))**2), &
                              & sqrt(1.0d0-(def_quantities_real_(34,2)/def_quantities_real_(33,2))**2)
  write(100,'(a24,1p,2e23.15)') 'sqrt(1-(Rz/Rx)^2)     = ', &
                              & sqrt(1.0d0-(def_quantities_real_(35,1)/def_quantities_real_(33,1))**2), &
                              & sqrt(1.0d0-(def_quantities_real_(35,2)/def_quantities_real_(33,2))**2)

  write(100,'(a24,1p,2e23.15)') 'proper_radius_x       = ', def_quantities_real_(36,1), def_quantities_real_(36,2)
  write(100,'(a24,1p,2e23.15)') 'proper_radius_y       = ', def_quantities_real_(37,1), def_quantities_real_(37,2)
  write(100,'(a24,1p,2e23.15)') 'proper_radius_z       = ', def_quantities_real_(38,1), def_quantities_real_(38,2)
  write(100,'(a24,1p,2e23.15)') 'Proper axis ratio y/x = ', def_quantities_real_(37,1)/def_quantities_real_(36,1),  &
                                                         &  def_quantities_real_(37,2)/def_quantities_real_(36,2)      
  write(100,'(a24,1p,2e23.15)') 'Proper axis ratio z/x = ', def_quantities_real_(38,1)/def_quantities_real_(36,1),  &
                                                         &  def_quantities_real_(38,2)/def_quantities_real_(36,2)
  write(100,'(a24,1p,2e23.15)') 'sqrt(1-(pRy/pRx)^2)   = ', &
                              & sqrt(1.0d0-(def_quantities_real_(37,1)/def_quantities_real_(36,1))**2), &
                              & sqrt(1.0d0-(def_quantities_real_(37,2)/def_quantities_real_(36,2))**2)
  write(100,'(a24,1p,2e23.15)') 'sqrt(1-(pRz/pRx)^2)   = ', &
                              & sqrt(1.0d0-(def_quantities_real_(38,1)/def_quantities_real_(36,1))**2), &
                              & sqrt(1.0d0-(def_quantities_real_(38,2)/def_quantities_real_(36,2))**2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'schwarz_radi_sph_km   = ', def_quantities_real_(39,1), def_quantities_real_(39,2)
  write(100,'(a24,1p,2e23.15)') 'coord_radius_x_km     = ', def_quantities_real_(40,1), def_quantities_real_(40,2)
  write(100,'(a24,1p,2e23.15)') 'coord_radius_y_km     = ', def_quantities_real_(41,1), def_quantities_real_(41,2)
  write(100,'(a24,1p,2e23.15)') 'coord_radius_z_km     = ', def_quantities_real_(42,1), def_quantities_real_(42,2)
  write(100,'(a24,1p,2e23.15)') 'proper_radius_x_km    = ', def_quantities_real_(43,1), def_quantities_real_(43,2)
  write(100,'(a24,1p,2e23.15)') 'proper_radius_y_km    = ', def_quantities_real_(44,1), def_quantities_real_(44,2)
  write(100,'(a24,1p,2e23.15)') 'proper_radius_z_km    = ', def_quantities_real_(45,1), def_quantities_real_(45,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'rho_c                 = ', def_quantities_real_(46,1), def_quantities_real_(46,2)
  write(100,'(a24,1p,2e23.15)') 'pre_c                 = ', def_quantities_real_(47,1), def_quantities_real_(47,2)
  write(100,'(a24,1p,2e23.15)') 'epsi_c                = ', def_quantities_real_(48,1), def_quantities_real_(48,2)
  write(100,'(a24,1p,2e23.15)') 'q_c                   = ', def_quantities_real_(49,1), def_quantities_real_(49,2)
  write(100,'(a24,1p,2e23.15)') 'rho_max               = ', def_quantities_real_(50,1), def_quantities_real_(50,2)
  write(100,'(a24,1p,2e23.15)') 'pre_max               = ', def_quantities_real_(51,1), def_quantities_real_(51,2)
  write(100,'(a24,1p,2e23.15)') 'epsi_max              = ', def_quantities_real_(52,1), def_quantities_real_(52,2)
  write(100,'(a24,1p,2e23.15)') 'q_max                 = ', def_quantities_real_(53,1), def_quantities_real_(53,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'rho_c_cgs             = ', def_quantities_real_(54,1), def_quantities_real_(54,2)
  write(100,'(a24,1p,2e23.15)') 'pre_c_cgs             = ', def_quantities_real_(55,1), def_quantities_real_(55,2)
  write(100,'(a24,1p,2e23.15)') 'epsi_c_cgs            = ', def_quantities_real_(56,1), def_quantities_real_(56,2)
  write(100,'(a24,1p,2e23.15)') 'q_c_cgs               = ', def_quantities_real_(57,1), def_quantities_real_(57,2)
  write(100,'(a24,1p,2e23.15)') 'rho_max_cgs           = ', def_quantities_real_(58,1), def_quantities_real_(58,2)
  write(100,'(a24,1p,2e23.15)') 'pre_max_cgs           = ', def_quantities_real_(59,1), def_quantities_real_(59,2)
  write(100,'(a24,1p,2e23.15)') 'epsi_max_cgs          = ', def_quantities_real_(60,1), def_quantities_real_(60,2)
  write(100,'(a24,1p,2e23.15)') 'q_max_cgs             = ', def_quantities_real_(61,1), def_quantities_real_(61,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'zrb_xp_plus           = ', def_quantities_real_(62,1), def_quantities_real_(62,2)
  write(100,'(a24,1p,2e23.15)') 'zrb_xp_minus          = ', def_quantities_real_(63,1), def_quantities_real_(63,2)
  write(100,'(a24,1p,2e23.15)') 'zrb_yp_plus           = ', def_quantities_real_(64,1), def_quantities_real_(64,2)
  write(100,'(a24,1p,2e23.15)') 'zrb_yp_minus          = ', def_quantities_real_(65,1), def_quantities_real_(65,2)
  write(100,'(a24,1p,2e23.15)') 'zrb_zp_plus           = ', def_quantities_real_(66,1), def_quantities_real_(66,2)
  write(100,'(a24,1p,2e23.15)') 'zrb_zp_minus          = ', def_quantities_real_(67,1), def_quantities_real_(67,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'dhdr_x                = ', def_quantities_real_(68,1), def_quantities_real_(68,2)
  write(100,'(a24,1p,2e23.15)') 'dhdr_y                = ', def_quantities_real_(69,1), def_quantities_real_(69,2)
  write(100,'(a24,1p,2e23.15)') 'dhdr_z                = ', def_quantities_real_(70,1), def_quantities_real_(70,2)
  write(100,'(a24,1p,2e23.15)') 'chi_cusp              = ', def_quantities_real_(71,1), def_quantities_real_(71,2)

  write(100,'(a72)') '========================================================================'
  write(100,'(a24,1p,2e23.15)') 'circ_shift_xy         = ', def_quantities_real_(80,1), def_quantities_real_(80,2)
  write(100,'(a24,1p,2e23.15)') 'circ_line_xy          = ', def_quantities_real_(72,1), def_quantities_real_(72,2)
  write(100,'(a24,1p,2e23.15)') 'circ_line_yz          = ', def_quantities_real_(73,1), def_quantities_real_(73,2)
  write(100,'(a24,1p,2e23.15)') 'circ_line_zx          = ', def_quantities_real_(74,1), def_quantities_real_(74,2)
  write(100,'(a24,1p,2e23.15)') 'circ_surf_xy          = ', def_quantities_real_(75,1), def_quantities_real_(75,2)
  write(100,'(a24,1p,2e23.15)') 'circ_surf_yz          = ', def_quantities_real_(76,1), def_quantities_real_(76,2)
  write(100,'(a24,1p,2e23.15)') 'circ_surf_zx          = ', def_quantities_real_(77,1), def_quantities_real_(77,2)
  write(100,'(a24,1p,2e23.15)') 'omespx                = ', def_matter_param_real_(13,1), def_matter_param_real_(13,2)
  write(100,'(a24,1p,2e23.15)') 'omespy                = ', def_matter_param_real_(14,1), def_matter_param_real_(14,2)
  write(100,'(a24,1p,2e23.15)') 'omespz                = ', def_matter_param_real_(15,1), def_matter_param_real_(15,2)

!  do ib  = 1, 3
!    do ia  = 1, 3
!      i=i+1   !70,75,80,85,90,95,100,105,110
!      write(100,'(a4,i1,a1,i1,a17,1p,2e23.15)') 'Iij(',ia,',',ib,')             = ', &
!      &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!    
!      i=i+1   !71,76,81,86,91,96,101,106,111
!      write(100,'(a4,i1,a1,i1,a17,1p,2e23.15)') 'Itf(',ia,',',ib,')             = ', &
!      &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!
!      i=i+1   !72,77,82,87,92,97,102,107,112
!      write(100,'(a7,i1,a1,i1,a14,1p,2e23.15)') 'dt1Itf(',ia,',',ib,')          = ', &
!      &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!
!      i=i+1   !73,78,83,88,93,98,103,108,113
!      write(100,'(a7,i1,a1,i1,a14,1p,2e23.15)') 'dt2Itf(',ia,',',ib,')          = ', &
!      &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!
!      i=i+1   !74,79,84,89,94,99,104,109,114
!      write(100,'(a7,i1,a1,i1,a14,1p,2e23.15)') 'dt3Itf(',ia,',',ib,')          = ', &
!      &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!    end do
!  end do
!
!  i=i+1   !115
!  write(100,'(a24,1p,2e23.15)')   'LGW                   = ', &
!  &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!
!  do ia  = 1, 3
!    i=i+1   !116,117,118
!    write(100,'(a5,i1,a18,1p,2e23.15)') 'dJdt(',ia,')               = ', &
!    &   def_quantities_real_(i,1), def_quantities_real_(i,2)
!  end do
!
!  i=i+1   !119
!  write(100,'(a24,1p,2e23.15)')   'hplus                 = ', &
!  &  def_quantities_real_(i,1), def_quantities_real_(i,2)
!
!  i=i+1   !120
!  write(100,'(a24,1p,2e23.15)')   'hcross                = ', &
!  &  def_quantities_real_(i,1), def_quantities_real_(i,2)

  close(100)
!
end subroutine printout_physq_BNS_all_mpt
