module def_kerr_schild
  use phys_constant, only : long
  implicit none
  real(long) :: reh_ks, rpop_ks, rpom_ks, rmbp_ks, rmbm_ks, rmsp_ks, rmsm_ks
  real(long) :: roerg_ks, AH_area_ks, irr_mass_ks, kerr_a(1:3)

  real(long), pointer  ::  alph_ks(:,:,:),     psi_ks(:,:,:)
  real(long), pointer  ::  bvxd_ks(:,:,:),     bvyd_ks(:,:,:),     bvzd_ks(:,:,:)
  real(long), pointer  ::  bvxu_ks(:,:,:),     bvyu_ks(:,:,:),     bvzu_ks(:,:,:)
  real(long), pointer  ::  hxxd_ks(:,:,:),     hxyd_ks(:,:,:),     hxzd_ks(:,:,:)
  real(long), pointer  ::  hyyd_ks(:,:,:),     hyzd_ks(:,:,:),     hzzd_ks(:,:,:)
!  real(long), pointer  ::  qxxd_ks(:,:,:),     qxyd_ks(:,:,:),     qxzd_ks(:,:,:)
!  real(long), pointer  ::  qyyd_ks(:,:,:),     qyzd_ks(:,:,:),     qzzd_ks(:,:,:)
!  real(long), pointer  ::  sou_ks(:,:,:,:)
!  real(long), pointer  ::  trk_ks(:,:,:)
!  real(long), pointer  ::  detg_ks(:,:,:), detg_grid_ks(:,:,:)
  real(long), pointer  ::  tfkij_grid_ks(:,:,:,:,:), tfkij_ks(:,:,:,:,:)
  real(long), pointer  ::  bcpsi_ks(:,:),      bcdpsi_ks(:,:,:)
  real(long), pointer  ::  bcalph_ks(:,:),     bcdalph_ks(:,:,:)
  real(long), pointer  ::  bcbvxd_ks(:,:),     bcdbvxd_ks(:,:,:)
  real(long), pointer  ::  bcbvyd_ks(:,:),     bcdbvyd_ks(:,:,:)
  real(long), pointer  ::  bcbvzd_ks(:,:),     bcdbvzd_ks(:,:,:)
  real(long), pointer  ::  bchxxd_ks(:,:),     bcdhxxd_ks(:,:,:)
  real(long), pointer  ::  bchxyd_ks(:,:),     bcdhxyd_ks(:,:,:)
  real(long), pointer  ::  bchxzd_ks(:,:),     bcdhxzd_ks(:,:,:)
  real(long), pointer  ::  bchyyd_ks(:,:),     bcdhyyd_ks(:,:,:)
  real(long), pointer  ::  bchyzd_ks(:,:),     bcdhyzd_ks(:,:,:)
  real(long), pointer  ::  bchzzd_ks(:,:),     bcdhzzd_ks(:,:,:)

  real(long), pointer  ::  obpsi_ks(:,:)
  real(long), pointer  ::  obalph_ks(:,:)
  real(long), pointer  ::  obbvxd_ks(:,:)
  real(long), pointer  ::  obbvyd_ks(:,:)
  real(long), pointer  ::  obbvzd_ks(:,:)
  real(long), pointer  ::  obhxxd_ks(:,:)
  real(long), pointer  ::  obhxyd_ks(:,:)
  real(long), pointer  ::  obhxzd_ks(:,:)
  real(long), pointer  ::  obhyyd_ks(:,:)
  real(long), pointer  ::  obhyzd_ks(:,:)
  real(long), pointer  ::  obhzzd_ks(:,:)

!  real(long), pointer  ::  Fxu_ks(:,:,:),      Fyu_ks(:,:,:),      Fzu_ks(:,:,:)
!  real(long), pointer  ::  Fxu_grid_ks(:,:,:), Fyu_grid_ks(:,:,:), Fzu_grid_ks(:,:,:)
end module def_kerr_schild
