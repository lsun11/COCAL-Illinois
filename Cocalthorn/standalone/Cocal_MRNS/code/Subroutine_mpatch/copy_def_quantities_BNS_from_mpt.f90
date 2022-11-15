subroutine copy_def_quantities_BNS_from_mpt(impt)
  use def_quantities
  use def_quantities_mpt
  implicit none
  integer :: i, impt, ia, ib
!
   admmass                = def_quantities_real_( 1,impt) 
   komarmass              = def_quantities_real_( 2,impt)
   komarmass_nc           = def_quantities_real_( 3,impt)
   restmass               = def_quantities_real_( 4,impt)
   propermass             = def_quantities_real_( 5,impt)
   angmom                 = def_quantities_real_( 6,impt)
   admmass_asymp          = def_quantities_real_( 7,impt)
   komarmass_asymp        = def_quantities_real_( 8,impt) 
   angmom_asymp           = def_quantities_real_( 9,impt)
   admmom_asymp(1)        = def_quantities_real_(10,impt)
   admmom_asymp(2)        = def_quantities_real_(11,impt) 
   admmom_asymp(3)        = def_quantities_real_(12,impt)
!
   T_kinene               = def_quantities_real_(13,impt) 
   W_gravene              = def_quantities_real_(14,impt) 
   P_intene               = def_quantities_real_(15,impt) 
   M_emfene               = def_quantities_real_(16,impt)
   M_torBene              = def_quantities_real_(17,impt)
   M_polBene              = def_quantities_real_(18,impt) 
   M_eleEene              = def_quantities_real_(19,impt)

   Virial                 = def_quantities_real_(20,impt)

   ToverW                 = def_quantities_real_(21,impt) 
   PoverW                 = def_quantities_real_(22,impt) 
   MoverW                 = def_quantities_real_(23,impt)
   MtorBoverW             = def_quantities_real_(24,impt) 
   MpolBoverW             = def_quantities_real_(25,impt) 
   MeleEoverW             = def_quantities_real_(26,impt) 

   I_inertia              = def_quantities_real_(27,impt) 

   gravmass_sph           = def_quantities_real_(28,impt) 
   restmass_sph           = def_quantities_real_(29,impt) 
   propermass_sph         = def_quantities_real_(30,impt) 
   MoverR_sph             = def_quantities_real_(31,impt) 

   schwarz_radi_sph       = def_quantities_real_(32,impt) 
   coord_radius_x         = def_quantities_real_(33,impt) 
   coord_radius_y         = def_quantities_real_(34,impt) 
   coord_radius_z         = def_quantities_real_(35,impt) 
   proper_radius_x        = def_quantities_real_(36,impt) 
   proper_radius_y        = def_quantities_real_(37,impt) 
   proper_radius_z        = def_quantities_real_(38,impt) 

   schwarz_radi_sph_km    = def_quantities_real_(39,impt) 
   coord_radius_x_km      = def_quantities_real_(40,impt) 
   coord_radius_y_km      = def_quantities_real_(41,impt) 
   coord_radius_z_km      = def_quantities_real_(42,impt) 
   proper_radius_x_km     = def_quantities_real_(43,impt) 
   proper_radius_y_km     = def_quantities_real_(44,impt) 
   proper_radius_z_km     = def_quantities_real_(45,impt)

   rho_c                  = def_quantities_real_(46,impt) 
   pre_c                  = def_quantities_real_(47,impt) 
   epsi_c                 = def_quantities_real_(48,impt) 
   q_c                    = def_quantities_real_(49,impt) 
   rho_max                = def_quantities_real_(50,impt) 
   pre_max                = def_quantities_real_(51,impt) 
   epsi_max               = def_quantities_real_(52,impt) 
   q_max                  = def_quantities_real_(53,impt) 

   rho_c_cgs              = def_quantities_real_(54,impt) 
   pre_c_cgs              = def_quantities_real_(55,impt) 
   epsi_c_cgs             = def_quantities_real_(56,impt)
   q_c_cgs                = def_quantities_real_(57,impt) 
   rho_max_cgs            = def_quantities_real_(58,impt) 
   pre_max_cgs            = def_quantities_real_(59,impt) 
   epsi_max_cgs           = def_quantities_real_(60,impt) 
   q_max_cgs              = def_quantities_real_(61,impt)

   zrb_xp_plus            = def_quantities_real_(62,impt) 
   zrb_xp_minus           = def_quantities_real_(63,impt) 
   zrb_yp_plus            = def_quantities_real_(64,impt) 
   zrb_yp_minus           = def_quantities_real_(65,impt) 
   zrb_zp_plus            = def_quantities_real_(66,impt) 
   zrb_zp_minus           = def_quantities_real_(67,impt) 

   dhdr_x                 = def_quantities_real_(68,impt) 
   dhdr_y                 = def_quantities_real_(69,impt) 
   dhdr_z                 = def_quantities_real_(70,impt) 
   chi_cusp               = def_quantities_real_(71,impt) 

   circ_line_xy           = def_quantities_real_(72,impt) 
   circ_line_yz           = def_quantities_real_(73,impt) 
   circ_line_zx           = def_quantities_real_(74,impt) 
   circ_surf_xy           = def_quantities_real_(75,impt) 
   circ_surf_yz           = def_quantities_real_(76,impt) 
   circ_surf_zx           = def_quantities_real_(77,impt) 

   ome_cgs                = def_quantities_real_(78,impt)
   qua_loc_spin           = def_quantities_real_(79,impt)
   circ_shift_xy          = def_quantities_real_(80,impt)
!
!   do ib = 1, 3
!     do ia = 1, 3
!       i=i+1; def_quantities_real_(i,impt) = Iij(ia, ib)
!       i=i+1; def_quantities_real_(i,impt) = Itf(ia, ib)
!       i=i+1; def_quantities_real_(i,impt) = dt1Itf(ia, ib)
!       i=i+1; def_quantities_real_(i,impt) = dt2Itf(ia, ib)
!       i=i+1; def_quantities_real_(i,impt) = dt3Itf(ia, ib)
!     end do
!   end do
!
!   def_quantities_real_(i,impt) = LGW
!
!   do ia = 1, 3
!     i=i+1; def_quantities_real_(i,impt) = dJdt(ia)
!   end do
!
!   def_quantities_real_(i,impt) = hplus
!   def_quantities_real_(i,impt) = hcross
!
!
end subroutine copy_def_quantities_BNS_from_mpt
