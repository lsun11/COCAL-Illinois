subroutine calc_physq_cgs_peos
  use phys_constant, only : g, c, solmas
  use def_matter_parameter, only : ome, radi
  use def_quantities, only : ome_cgs, rho_c, pre_c, epsi_c, q_c, & 
  &          rho_c_cgs, pre_c_cgs, epsi_c_cgs, q_c_cgs, &
  &          rho_max, pre_max, epsi_max, q_max,   &
  &          rho_max_cgs, pre_max_cgs, epsi_max_cgs, q_max_cgs,   &
  &          schwarz_radi_sph, schwarz_radi_sph_km, &
  &          coord_radius_x, coord_radius_y, coord_radius_z, &
  &          proper_radius_x, proper_radius_y, proper_radius_z, &
  &          coord_radius_x_km, coord_radius_y_km, coord_radius_z_km, &
  &          proper_radius_x_km, proper_radius_y_km, proper_radius_z_km
  implicit none
  real(8) :: MM = solmas, LL = g*solmas/c**2, TT = g*solmas/c**3
!
!  rhocgs_c = rho_c*c**6/(g**3*solmas**2)
!  precgs_c = pre_c*c**8/(g**3*solmas**2)
!  epsicgs_c = epsi_c*c**6/(g**3*solmas**2)
!  radicgs = radi*g*solmas*1.0d-5/c**2
!
  rho_c_cgs = rho_c*MM/LL**3
  pre_c_cgs = pre_c*MM/(LL*TT**2)
  epsi_c_cgs = epsi_c*MM/LL**3
  q_c_cgs = pre_c_cgs/rho_c_cgs
!
  rho_max_cgs = rho_max*MM/LL**3
  pre_max_cgs = pre_max*MM/(LL*TT**2)
  epsi_max_cgs = epsi_max*MM/LL**3
  q_max_cgs = pre_max_cgs/rho_max_cgs
!
  schwarz_radi_sph_km = schwarz_radi_sph*LL*1.0d-5
  coord_radius_x_km = coord_radius_x*LL*1.0d-5
  coord_radius_y_km = coord_radius_y*LL*1.0d-5
  coord_radius_z_km = coord_radius_z*LL*1.0d-5
  proper_radius_x_km = proper_radius_x*LL*1.0d-5
  proper_radius_y_km = proper_radius_y*LL*1.0d-5
  proper_radius_z_km = proper_radius_z*LL*1.0d-5
!
  ome_cgs = ome/radi/TT
end subroutine calc_physq_cgs_peos
