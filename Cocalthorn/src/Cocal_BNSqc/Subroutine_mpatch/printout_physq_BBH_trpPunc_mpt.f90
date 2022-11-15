subroutine printout_physq_BBH_trpPunc_mpt(iseq,impt)
  use phys_constant, only : long, nmpt
  use grid_parameter, only : rgin
  use def_quantities
  use def_bh_parameter
  use def_binary_parameter, only : sepa, dis
  use def_quantities_bh
  implicit none
  real(long) :: fixedvir
  integer, intent(in) :: iseq, impt
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/)
!
  if (iseq.eq.1) then 
    open(110,file='bbhphyseq_mpt'//np(impt)//'.dat',status='unknown')
  end if
  write(110,*) '== Sequence number == ', iseq
  write(110,*) '## BH Radii in G = c = 1 unit ##'
  write(110,'(a22,1p,2e23.15)') ' Radius & dis of  BH= ', rgin, dis
  write(110,'(a22,1p,2e23.15)') ' Separation         = ', sepa
!
  write(110,*) '## BH orbital and spin angular velocity parameters ##'
  write(110,'(a22,1p,2e23.15)') ' Omega              = ', ome_bh
  write(110,'(a22,1p,2e23.15)') ' Spin      of  BH   = ', spin_bh
  write(110,'(a22,1p,2e23.15)') ' Spin axis of  BH   = ', th_spin_bh_deg, &
  &                                                       phi_spin_bh_deg
  write(110,'(a22,1p,2e23.15)') ' Puncture spin Sx   = ', spin_pBH(1)
  write(110,'(a22,1p,2e23.15)') ' Puncture spin Sy   = ', spin_pBH(2)
  write(110,'(a22,1p,2e23.15)') ' Puncture spin Sz   = ', spin_pBH(3)
  write(110,*) '## BH linear momentum parameters ##'
  write(110,'(a22,1p,2e23.15)') ' Px asympto & punct = ', admmom_asymp(1), &
  &                                                       mom_pBH(1)
  write(110,'(a22,1p,2e23.15)') ' Py asympto & punct = ', admmom_asymp(2), &
  &                                                       mom_pBH(2)
  write(110,'(a22,1p,2e23.15)') ' Pz asympto & punct = ', admmom_asymp(3), &
  &                                                       mom_pBH(3)
!
  write(110,*) '## Mass and angular momentum (G = c = 1 unit) ##'
  write(110,'(a22,1p,2e23.15)') ' Puncture mass para = ', mass_pBH
  write(110,'(a22,1p,2e23.15)') ' M_ADM (asymp & vol)= ', admmass_asymp, &
  &                                                       admmass
  write(110,'(a22,1p,2e23.15)') ' M_K   (asymp & vol)= ', komarmass_asymp, &
  &                                                       komarmass
  write(110,'(a22,1p,2e23.15)') ' J     (asymptotic) = ', angmom_asymp
  write(110,'(a22,1p,2e23.15)') ' AH mass of BH      = ', AHmass
  write(110,'(a22,1p,2e23.15)') ' AH area of BH      = ', AHarea
  write(110,'(a22,1p,2e23.15)') ' AH spin of BH      = ', AHspin
!
  fixedvir  = (admmass_asymp - komarmass_asymp)/admmass_asymp
  write(110,'(a22,1p,2e23.15)') ' 1 - M_K/M_ADM      = ', fixedvir
!
end subroutine printout_physq_BBH_trpPunc_mpt
