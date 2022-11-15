subroutine calc_Mirr_Mchr_TW
  use phys_constant, only  : long, pi
  use def_bht_parameter, only : am_bh
  use def_quantities, only : admmass_asymp, app_hor_area_bh, irredmass, &
       &                     bindingene, christmass, angmom_asymp, propermass, &
       &                     W_gravene_omeJ, T_kinene_omeJ, ToverW_omeJ
  implicit none
  real(long) :: spinpar
!
  irredmass = dsqrt(app_hor_area_bh/(16.0d0*pi))    

  bindingene = admmass_asymp - irredmass

  spinpar = am_bh/(irredmass*irredmass)
  christmass = irredmass*dsqrt(1.0d0 + spinpar*spinpar/4.0d0)
!
  W_gravene_omeJ = admmass_asymp - christmass - propermass - T_kinene_omeJ
  ToverW_omeJ = T_kinene_omeJ/dabs(W_gravene_omeJ)

end subroutine calc_Mirr_Mchr_TW
