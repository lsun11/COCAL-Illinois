subroutine printout_physq_BBH_mpt(iseq)
  use phys_constant, only : long, nmpt
  use grid_parameter, only : rgin
  use def_quantities
  use def_bh_parameter
  use def_binary_parameter, only : sepa, dis
  use def_quantities_bh
  implicit none
  real(long) :: fixedvir
  integer, intent(in) :: iseq
  integer :: impt
!
  if (iseq.eq.1) then 
    open(110,file='bbhphyseq.dat',status='unknown')
  end if
  write(110,*) '== Sequence number == ', iseq
  write(110,*) '## BH Radii in G = c = 1 unit ##'
  call copy_grid_parameter_from_mpt(1)
  call copy_def_binary_parameter_from_mpt(1)
  write(110,'(a22,1p,2e23.15)') ' Radius & dis 1st BH= ', rgin, dis
  call copy_grid_parameter_from_mpt(2)
  call copy_def_binary_parameter_from_mpt(2)
  write(110,'(a22,1p,2e23.15)') ' Radius & dis 2nd BH= ', rgin, dis
  call copy_def_binary_parameter_from_mpt(1)
  write(110,'(a22,1p,2e23.15)') ' Separation         = ', sepa
!
  write(110,*) '## BBH orbital and spin angular velocity parameters ##'
  call copy_def_bh_parameter_from_mpt(1)
  write(110,'(a22,1p,2e23.15)') ' Omega              = ', ome_bh
  write(110,'(a22,1p,2e23.15)') ' Spin      1st BH   = ', spin_bh
  write(110,'(a22,1p,2e23.15)') ' Spin axis 1st BH   = ', th_spin_bh_deg, &
  &                                                       phi_spin_bh_deg
  call copy_def_bh_parameter_from_mpt(2)
  write(110,'(a22,1p,2e23.15)') ' Spin      2nd BH   = ', spin_bh
  write(110,'(a22,1p,2e23.15)') ' Spin axis 2nd BH   = ', th_spin_bh_deg, &
  &                                                       phi_spin_bh_deg
  write(110,*) '## Mass and angular momentum (G = c = 1 unit) ##'
  call copy_def_quantities_from_mpt(nmpt)
  write(110,'(a22,1p,2e23.15)') ' M_ADM (asymptotic) = ', admmass_asymp
  write(110,'(a22,1p,2e23.15)') ' M_K Komar mass     = ', komarmass_asymp
  write(110,'(a22,1p,2e23.15)') ' J     (asymptotic) = ', angmom_asymp
  call copy_def_quantities_bh_from_mpt(1)
  write(110,'(a22,1p,2e23.15)') ' AH mass   1st BH   = ', AHmass
  write(110,'(a22,1p,2e23.15)') ' AH area   1st BH   = ', AHarea
  write(110,'(a22,1p,2e23.15)') ' AH spin   1st BH   = ', AHspin
  call copy_def_quantities_bh_from_mpt(2)
  write(110,'(a22,1p,2e23.15)') ' AH mass   2nd BH   = ', AHmass
  write(110,'(a22,1p,2e23.15)') ' AH area   2nd BH   = ', AHarea
  write(110,'(a22,1p,2e23.15)') ' AH spin   2nd BH   = ', AHspin
!
  fixedvir  = (admmass_asymp - komarmass_asymp)/admmass_asymp
  write(110,'(a22,1p,2e23.15)') ' 1 - M_K/M_ADM      = ', fixedvir
!
end subroutine printout_physq_BBH_mpt
