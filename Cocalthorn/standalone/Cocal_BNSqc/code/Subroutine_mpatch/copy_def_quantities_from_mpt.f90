subroutine copy_def_quantities_from_mpt(impt)
  use def_quantities
  use def_quantities_mpt
  implicit none
  integer :: i, impt, ia, ib
!
  i=0
   i=i+1; admmass                = def_quantities_real_(i,impt)  !1
   i=i+1; komarmass              = def_quantities_real_(i,impt)  !2
   i=i+1; komarmass_nc           = def_quantities_real_(i,impt)  !3
   i=i+1; restmass               = def_quantities_real_(i,impt)  !4
   i=i+1; propermass             = def_quantities_real_(i,impt)  !5
   i=i+1; angmom                 = def_quantities_real_(i,impt)  !6
   i=i+1; admmass_asymp          = def_quantities_real_(i,impt)  !7
   i=i+1; komarmass_asymp        = def_quantities_real_(i,impt)  !8
   i=i+1; angmom_asymp           = def_quantities_real_(i,impt)  !9
!
   i=i+1; admmom_asymp(1)        = def_quantities_real_(i,impt)  !10
   i=i+1; admmom_asymp(2)        = def_quantities_real_(i,impt)  !11
   i=i+1; admmom_asymp(3)        = def_quantities_real_(i,impt)  !12
!
   i=i+1; T_kinene               = def_quantities_real_(i,impt)  !13
   i=i+1; W_gravene              = def_quantities_real_(i,impt)  !14
   i=i+1; P_intene               = def_quantities_real_(i,impt)  !15
   i=i+1; M_emfene               = def_quantities_real_(i,impt)  !16
   i=i+1; M_torBene              = def_quantities_real_(i,impt)  !17
   i=i+1; M_polBene              = def_quantities_real_(i,impt)  !18
   i=i+1; M_eleEene              = def_quantities_real_(i,impt)  !19
   i=i+1; Virial                 = def_quantities_real_(i,impt)  !20
   i=i+1; ToverW                 = def_quantities_real_(i,impt)  !21
   i=i+1; PoverW                 = def_quantities_real_(i,impt)  !22
   i=i+1; MoverW                 = def_quantities_real_(i,impt)  !23
   i=i+1; MtorBoverW             = def_quantities_real_(i,impt)  !24
   i=i+1; MpolBoverW             = def_quantities_real_(i,impt)  !25
   i=i+1; MeleEoverW             = def_quantities_real_(i,impt)  !26
   i=i+1; I_inertia              = def_quantities_real_(i,impt)  !27
   i=i+1; gravmass_sph           = def_quantities_real_(i,impt)  !28
   i=i+1; restmass_sph           = def_quantities_real_(i,impt)  !29
   i=i+1; propermass_sph         = def_quantities_real_(i,impt)  !30
   i=i+1; MoverR_sph             = def_quantities_real_(i,impt)  !31
   i=i+1; schwarz_radi_sph       = def_quantities_real_(i,impt)  !32
   i=i+1; schwarz_radi_sph_km    = def_quantities_real_(i,impt)  !33
   i=i+1; coord_radius_x         = def_quantities_real_(i,impt)  !34
   i=i+1; coord_radius_y         = def_quantities_real_(i,impt)  !35
   i=i+1; coord_radius_z         = def_quantities_real_(i,impt)  !36
   i=i+1; proper_radius_x        = def_quantities_real_(i,impt)  !37
   i=i+1; proper_radius_y        = def_quantities_real_(i,impt)  !38
   i=i+1; proper_radius_z        = def_quantities_real_(i,impt)  !39
   i=i+1; rho_c                  = def_quantities_real_(i,impt)  !40
   i=i+1; pre_c                  = def_quantities_real_(i,impt)  !41
   i=i+1; q_c                    = def_quantities_real_(i,impt)  !42
   i=i+1; rho_max                = def_quantities_real_(i,impt)  !43
   i=i+1; pre_max                = def_quantities_real_(i,impt)  !44
   i=i+1; epsi_max               = def_quantities_real_(i,impt)  !45
   i=i+1; q_max                  = def_quantities_real_(i,impt)  !46
   i=i+1; coord_radius_x_km      = def_quantities_real_(i,impt)  !47
   i=i+1; coord_radius_y_km      = def_quantities_real_(i,impt)  !48
   i=i+1; coord_radius_z_km      = def_quantities_real_(i,impt)  !49
   i=i+1; proper_radius_x_km     = def_quantities_real_(i,impt)  !50
   i=i+1; proper_radius_y_km     = def_quantities_real_(i,impt)  !51
   i=i+1; proper_radius_z_km     = def_quantities_real_(i,impt)  !52
   i=i+1; rho_c_cgs              = def_quantities_real_(i,impt)  !53
   i=i+1; pre_c_cgs              = def_quantities_real_(i,impt)  !54
   i=i+1; epsi_c_cgs             = def_quantities_real_(i,impt)  !55
   i=i+1; q_c_cgs                = def_quantities_real_(i,impt)  !56
   i=i+1; rho_max_cgs            = def_quantities_real_(i,impt)  !57
   i=i+1; pre_max_cgs            = def_quantities_real_(i,impt)  !58
   i=i+1; epsi_max_cgs           = def_quantities_real_(i,impt)  !59
   i=i+1; q_max_cgs              = def_quantities_real_(i,impt)  !60
   i=i+1; zrb_xp_plus            = def_quantities_real_(i,impt)  !61
   i=i+1; zrb_xp_minus           = def_quantities_real_(i,impt)  !62
   i=i+1; zrb_yp_plus            = def_quantities_real_(i,impt)  !63
   i=i+1; zrb_yp_minus           = def_quantities_real_(i,impt)  !64
   i=i+1; zrb_zp_plus            = def_quantities_real_(i,impt)  !65
   i=i+1; zrb_zp_minus           = def_quantities_real_(i,impt)  !66
   i=i+1; dhdr_x                 = def_quantities_real_(i,impt)  !67
   i=i+1; dhdr_y                 = def_quantities_real_(i,impt)  !68
   i=i+1; dhdr_z                 = def_quantities_real_(i,impt)  !69
!
   do ib = 1, 3
     do ia = 1, 3
       i=i+1; Iij(ia, ib)        = def_quantities_real_(i,impt)  !70,75,80,85,90,95,100,105,110
       i=i+1; Itf(ia, ib)        = def_quantities_real_(i,impt)  !71,76,81,86,91,96,101,106,111
       i=i+1; dt1Itf(ia, ib)     = def_quantities_real_(i,impt)  !72,77,82,87,92,97,102,107,112
       i=i+1; dt2Itf(ia, ib)     = def_quantities_real_(i,impt)  !73,78,83,88,93,98,103,108,113
       i=i+1; dt3Itf(ia, ib)     = def_quantities_real_(i,impt)  !74,79,84,89,94,99,104,109,114
     end do
   end do
!
   i=i+1; LGW                    = def_quantities_real_(i,impt)  !115
!
   do ia = 1, 3
     i=i+1; dJdt(ia)             = def_quantities_real_(i,impt)  !116,117,118
   end do
!
   i=i+1; hplus                  = def_quantities_real_(i,impt)  !119
   i=i+1; hcross                 = def_quantities_real_(i,impt)  !120
!
   i=i+1; chi_cusp               = def_quantities_real_(i,impt)  !121
!
   i=i+1; ome_cgs                = def_quantities_real_(i,impt)  !122
   i=i+1; qua_loc_spin           = def_quantities_real_(i,impt)  !123
!
end subroutine copy_def_quantities_from_mpt
