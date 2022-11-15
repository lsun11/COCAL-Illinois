subroutine copy_def_quantities_to_mpt(impt)
  use def_quantities
  use def_quantities_mpt
  implicit none
  integer :: i, impt, ia, ib
!
  i=0
   i=i+1; def_quantities_real_(i,impt) = admmass
   i=i+1; def_quantities_real_(i,impt) = komarmass
   i=i+1; def_quantities_real_(i,impt) = komarmass_nc
   i=i+1; def_quantities_real_(i,impt) = restmass
   i=i+1; def_quantities_real_(i,impt) = propermass
   i=i+1; def_quantities_real_(i,impt) = angmom
   i=i+1; def_quantities_real_(i,impt) = admmass_asymp
   i=i+1; def_quantities_real_(i,impt) = komarmass_asymp
   i=i+1; def_quantities_real_(i,impt) = angmom_asymp
!
   i=i+1; def_quantities_real_(i,impt) = admmom_asymp(1)
   i=i+1; def_quantities_real_(i,impt) = admmom_asymp(2)
   i=i+1; def_quantities_real_(i,impt) = admmom_asymp(3)
!
   i=i+1; def_quantities_real_(i,impt) = T_kinene
   i=i+1; def_quantities_real_(i,impt) = W_gravene
   i=i+1; def_quantities_real_(i,impt) = P_intene
   i=i+1; def_quantities_real_(i,impt) = M_emfene
   i=i+1; def_quantities_real_(i,impt) = M_torBene
   i=i+1; def_quantities_real_(i,impt) = M_polBene
   i=i+1; def_quantities_real_(i,impt) = M_eleEene
   i=i+1; def_quantities_real_(i,impt) = Virial
   i=i+1; def_quantities_real_(i,impt) = ToverW
   i=i+1; def_quantities_real_(i,impt) = PoverW
   i=i+1; def_quantities_real_(i,impt) = MoverW
   i=i+1; def_quantities_real_(i,impt) = MtorBoverW
   i=i+1; def_quantities_real_(i,impt) = MpolBoverW
   i=i+1; def_quantities_real_(i,impt) = MeleEoverW
   i=i+1; def_quantities_real_(i,impt) = I_inertia
   i=i+1; def_quantities_real_(i,impt) = gravmass_sph
   i=i+1; def_quantities_real_(i,impt) = restmass_sph
   i=i+1; def_quantities_real_(i,impt) = propermass_sph
   i=i+1; def_quantities_real_(i,impt) = MoverR_sph
   i=i+1; def_quantities_real_(i,impt) = schwarz_radi_sph
   i=i+1; def_quantities_real_(i,impt) = schwarz_radi_sph_km
   i=i+1; def_quantities_real_(i,impt) = coord_radius_x
   i=i+1; def_quantities_real_(i,impt) = coord_radius_y
   i=i+1; def_quantities_real_(i,impt) = coord_radius_z
   i=i+1; def_quantities_real_(i,impt) = proper_radius_x
   i=i+1; def_quantities_real_(i,impt) = proper_radius_y
   i=i+1; def_quantities_real_(i,impt) = proper_radius_z
   i=i+1; def_quantities_real_(i,impt) = rho_c
   i=i+1; def_quantities_real_(i,impt) = pre_c
   i=i+1; def_quantities_real_(i,impt) = q_c
   i=i+1; def_quantities_real_(i,impt) = rho_max
   i=i+1; def_quantities_real_(i,impt) = pre_max
   i=i+1; def_quantities_real_(i,impt) = epsi_max
   i=i+1; def_quantities_real_(i,impt) = q_max
   i=i+1; def_quantities_real_(i,impt) = coord_radius_x_km
   i=i+1; def_quantities_real_(i,impt) = coord_radius_y_km
   i=i+1; def_quantities_real_(i,impt) = coord_radius_z_km
   i=i+1; def_quantities_real_(i,impt) = proper_radius_x_km
   i=i+1; def_quantities_real_(i,impt) = proper_radius_y_km
   i=i+1; def_quantities_real_(i,impt) = proper_radius_z_km
   i=i+1; def_quantities_real_(i,impt) = rho_c_cgs
   i=i+1; def_quantities_real_(i,impt) = pre_c_cgs
   i=i+1; def_quantities_real_(i,impt) = epsi_c_cgs
   i=i+1; def_quantities_real_(i,impt) = q_c_cgs
   i=i+1; def_quantities_real_(i,impt) = rho_max_cgs
   i=i+1; def_quantities_real_(i,impt) = pre_max_cgs
   i=i+1; def_quantities_real_(i,impt) = epsi_max_cgs
   i=i+1; def_quantities_real_(i,impt) = q_max_cgs
   i=i+1; def_quantities_real_(i,impt) = zrb_xp_plus
   i=i+1; def_quantities_real_(i,impt) = zrb_xp_minus
   i=i+1; def_quantities_real_(i,impt) = zrb_yp_plus
   i=i+1; def_quantities_real_(i,impt) = zrb_yp_minus
   i=i+1; def_quantities_real_(i,impt) = zrb_zp_plus
   i=i+1; def_quantities_real_(i,impt) = zrb_zp_minus
   i=i+1; def_quantities_real_(i,impt) = dhdr_x
   i=i+1; def_quantities_real_(i,impt) = dhdr_y
   i=i+1; def_quantities_real_(i,impt) = dhdr_z
!
   do ib = 1, 3
     do ia = 1, 3
       i=i+1; def_quantities_real_(i,impt) = Iij(ia, ib)
       i=i+1; def_quantities_real_(i,impt) = Itf(ia, ib)
       i=i+1; def_quantities_real_(i,impt) = dt1Itf(ia, ib)
       i=i+1; def_quantities_real_(i,impt) = dt2Itf(ia, ib)
       i=i+1; def_quantities_real_(i,impt) = dt3Itf(ia, ib)
     end do
   end do
!
   i=i+1; def_quantities_real_(i,impt) = LGW
!
   do ia = 1, 3
     i=i+1; def_quantities_real_(i,impt) = dJdt(ia)
   end do
!
   i=i+1; def_quantities_real_(i,impt) = hplus
   i=i+1; def_quantities_real_(i,impt) = hcross
!
   i=i+1; def_quantities_real_(i,impt) = chi_cusp
!
   i=i+1; def_quantities_real_(i,impt) = ome_cgs
   i=i+1; def_quantities_real_(i,impt) = qua_loc_spin
!
end subroutine copy_def_quantities_to_mpt
