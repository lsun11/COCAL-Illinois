subroutine printout_physq(iseq,flag_restmass)
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, ntfpolp, sw_mass_iter
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome, radi, pinx
  use def_quantities
  implicit none
  real(8) :: fixeddlm, fixedvir
  integer, intent(in) :: iseq, flag_restmass
!
  if (iseq.eq.1) then 
    open(100,file='rnsphyseq.dat',status='unknown')
  end if
  if (flag_restmass.le.1) write(100,*) '#### Solution did not converge #### '
  write(100,*) '== Sequence number == ', iseq
  write(100,*) '## Coordinate and Proper Radii in K = 1 unit ##'
  write(100,'(a22,1p,2e23.15)') ' NS radius along x  = ', &
  &                              coord_radius_x, proper_radius_x
  write(100,'(a22,1p,2e23.15)') ' NS radius along y  = ', &
  &                              coord_radius_y, proper_radius_y
  write(100,'(a22,1p,2e23.15)') ' NS radius along z  = ', &
  &                              coord_radius_z, proper_radius_z
  write(100,'(a22,1p,2e23.15)') ' Axis ratio y/x     = ', &
  &           coord_radius_y/coord_radius_x, proper_radius_y/proper_radius_x
  write(100,'(a22,1p,2e23.15)') ' Axis ratio z/x     = ', &
  &           coord_radius_z/coord_radius_x, proper_radius_z/proper_radius_x
  write(100,*) '## Radii in Req = 1 unit ##'
  write(100,'(a22,1p,2e23.15)') ' NS radius along x  = ', rs(ntfeq,npfxzp)
  write(100,'(a22,1p,2e23.15)') ' NS radius along y  = ', rs(ntfeq,npfyzp)
  write(100,'(a22,1p,2e23.15)') ' NS radius along z  = ', rs(ntfpolp,npfxzp)
  write(100,*) '## Eccentricity in coordinate and Proper Radii ##'
  write(100,'(a22,1p,2e23.15)') ' sqrt(1-(Rz/Rx)^2)  = ', &
  &           sqrt(1.0d0-(coord_radius_z/coord_radius_x)**2), &
  &           sqrt(1.0d0-(proper_radius_z/proper_radius_x)**2)
!
  write(100,*) '## M = spherical gravitational mass (K = 1 unit) ##'
  write(100,'(a22,1p,2e23.15)') ' Omega M and Omega  = ', & 
  &                               ome/radi*gravmass_sph, ome/radi
  write(100,'(a22,1p,2e23.15)') ' M_ADM              = ', admmass
  write(100,'(a22,1p,2e23.15)') ' M_K Komar mass     = ', komarmass
  write(100,'(a22,1p,2e23.15)') ' Rest mass M_0 1star= ', restmass
  write(100,'(a22,1p,2e23.15)') ' Proper mass M_p    = ', propermass
  write(100,'(a22,1p,2e23.15)') ' Angular momentum J = ', angmom, angmom_asymp
  write(100,'(a22,1p,2e23.15)') ' Spherical rest mass= ', restmass_sph
  write(100,'(a22,1p,2e23.15)') ' Spherical grav mass= ', gravmass_sph
  write(100,'(a22,1p,2e23.15)') ' Spherical M/R      = ', MoverR_sph
  write(100,'(a22,1p,2e23.15)') ' Schwarzschildradius= ', schwarz_radi_sph
!
  write(100,'(a22,1p,2e23.15)') ' E/M = (M_ADM-M)/M  = ', &
  &                                            admmass/gravmass_sph - 1.0d0
  write(100,'(a22,1p,2e23.15)') ' J/M^2 and J/M_ADM^2= ', &
  &                     angmom/gravmass_sph**2.0d0, angmom/admmass**2.0d0
!
  write(100,'(a22,1p,2e23.15)') ' T kinetic energy   = ', T_kinene
  write(100,'(a22,1p,2e23.15)') ' W grav energy      = ', W_gravene
  write(100,'(a22,1p,2e23.15)') ' T/|W|              = ', ToverW
  write(100,'(a22,1p,2e23.15)') ' Moment of inertia  = ', I_inertia
!
  write(100,*) '## Virial and Rest mass accuracy ##'
  fixedvir  = (admmass  -    komarmass)/admmass
  fixeddlm  = (restmass - restmass_sph)/restmass_sph
  write(100,'(a22,1p,2e23.15)') ' 1 - M_K/M_ADM      = ', fixedvir
  write(100,'(a22,1p,2e23.15)') ' M_0/M_0sph - 1     = ', fixeddlm
!
  write(100,*) '## Maximums of the density, pressure and p/rho ##'
  write(100,'(a22,1p,2e23.15)') ' rho_max            = ', rho_max
  write(100,'(a22,1p,2e23.15)') ' pre_max            = ', pre_max
  write(100,'(a22,1p,2e23.15)') ' epsilon_max        = ', epsi_max
  write(100,'(a22,1p,2e23.15)') ' (p/rho)_max        = ', q_max
!
  write(100,*) '## Red and blue shift ##'
  write(100,'(a22,1p,2e23.15)') ' surface on x axis  = ', &
  &                               zrb_xp_plus, zrb_xp_minus 
  write(100,'(a22,1p,2e23.15)') ' surface on y axis  = ', &
  &                               zrb_yp_plus, zrb_yp_minus
  write(100,'(a22,1p,2e23.15)') ' surface on z axis  = ', &
  &                               zrb_zp_plus, zrb_zp_minus
!
end subroutine printout_physq
