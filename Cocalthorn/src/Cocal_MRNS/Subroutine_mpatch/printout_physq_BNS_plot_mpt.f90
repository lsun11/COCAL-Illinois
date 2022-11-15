subroutine printout_physq_BNS_plot_mpt(filename,impt)
  use phys_constant,  only : nmpt, g, c, solmas, long
  use grid_parameter, only : ntfeq, npfxzp, npfyzp, ntfpolp, sw_mass_iter
  use def_matter, only : emd, rs
  use def_matter_parameter, only : ome, radi, pinx
  use def_quantities
  use def_quantities_mpt
  use def_matter_parameter_mpt
  use grid_parameter_mpt
  use def_binary_parameter
  use def_binary_parameter_mpt
  use make_array_1d
  implicit none
  real(long), pointer :: printmat(:) 
  real(8) :: MM = solmas, LL = g*solmas/c**2, TT = g*solmas/c**3
  real(8) :: Msph_tot
  real(8) :: fixeddlm, fixedvir
  integer :: i,ia,ib
  integer, intent(in) :: impt
  character(len=30), intent(in) :: filename
!
  if ((impt.ne.1).and.(impt.ne.2))  then
    write(6,*) "printout_physq_BNS_plot_mpt: Need patch 1 or 2...exiting"
    write(6,*) "impt=",impt
    stop
  end if

  Msph_tot = def_quantities_real_(28,1) + def_quantities_real_(28,2)
  call alloc_array1d(printmat, 1, 100)

  printmat( 1) = def_matter_param_real_(3,impt)    ! Omega R0  !
  printmat( 2) = def_matter_param_real_(5,impt)    ! R0     !
  printmat( 3) = def_matter_param_real_(3,impt)/def_matter_param_real_(5,impt) !Omega!
  printmat( 4) = def_quantities_real_(122,impt)    ! Omega [rad/sec] !
  printmat( 5) = def_matter_param_real_(3,impt)/def_matter_param_real_(5,impt)*Msph_tot   ! Omega M  !
  printmat( 6) = sepa_(impt)*def_matter_param_real_(5,impt)            ! Distance d_s =sepa*R0 !
  printmat( 7) = sepa_(impt)*def_matter_param_real_(5,impt)*LL*1.0d-5  ! Distance d_s [km]!
  printmat( 8) = dis_(impt)*def_matter_param_real_(5,impt)             ! Distance d =dis*R0!
  printmat( 9) = dis_(impt)*def_matter_param_real_(5,impt)*LL*1.0d-5   ! Distance d [km]!
  printmat(10) = def_quantities_real_(4,impt)     ! restmass  !
  printmat(11) = def_quantities_real_(5,impt)     ! propermass  !
  printmat(12) = def_quantities_real_(6,impt)     ! angmom   !
  printmat(13) = (def_quantities_real_(4,impt) - def_quantities_real_(29,impt))/def_quantities_real_(29,impt)   ! M_0/M_0sph - 1   !
  printmat(14) = def_quantities_real_(7,nmpt)  ! admmass_asymp  !
  printmat(15) = def_quantities_real_(8,nmpt)  ! komarmass_asymp!
  printmat(16) = def_quantities_real_(9,nmpt)  ! angmom_asymp     !
  printmat(17) = def_quantities_real_(10,nmpt) ! admmom_asymp(1)  !
  printmat(18) = def_quantities_real_(11,nmpt) ! admmom_asymp(2)  !
  printmat(19) = def_quantities_real_(12,nmpt) ! admmom_asymp(3)  !
  printmat(20) = def_quantities_real_(7,nmpt)/Msph_tot - 1.0d0    ! E/M = (M_ADM-M)/M   !
  printmat(21) = def_quantities_real_(9,nmpt)/Msph_tot**2.0d0     ! J/M^2  !
  printmat(22) = def_quantities_real_(9,nmpt)/def_quantities_real_(7,nmpt)**2.0d0   ! J/M_ADM^2  !
  printmat(23) = (def_quantities_real_(7,nmpt) - def_quantities_real_(8,nmpt))/def_quantities_real_(7,nmpt)   ! 1 - M_K/M_ADM !
  printmat(24) = def_quantities_real_(28,impt)     ! gravmass_sph !
  printmat(25) = def_quantities_real_(29,impt)     ! restmass_sph !
  printmat(26) = def_quantities_real_(30,impt)     ! propermass_sph !
  printmat(27) = def_quantities_real_(31,impt)     ! MoverR_sph  !
  printmat(28) = def_quantities_real_(32,impt)     ! schwarz_radi_sph  !
  printmat(29) = def_quantities_real_(33,impt)     ! schwarz_radi_sph_km !
  printmat(30) = def_quantities_real_(34,impt)     ! coord_radius_x  !
  printmat(31) = def_quantities_real_(35,impt)     ! coord_radius_y  !
  printmat(32) = def_quantities_real_(36,impt)     ! coord_radius_z  !
  printmat(33) = def_quantities_real_(35,impt)/def_quantities_real_(34,impt)   ! Coord axis ratio y/x   !
  printmat(34) = def_quantities_real_(36,impt)/def_quantities_real_(34,impt)   ! Coord axis ratio z/x  !
  printmat(35) = sqrt(1.0d0-(def_quantities_real_(36,impt)/def_quantities_real_(34,impt))**2)  !sqrt(1-(Rz/Rx)^2)  !
  printmat(36) = def_quantities_real_(37,impt)     ! proper_radius_x !
  printmat(37) = def_quantities_real_(38,impt)     ! proper_radius_y !
  printmat(38) = def_quantities_real_(39,impt)     ! proper_radius_z !
  printmat(39) = def_quantities_real_(38,impt)/def_quantities_real_(37,impt)  ! Proper axis ratio y/x  !
  printmat(40) = def_quantities_real_(39,impt)/def_quantities_real_(37,impt)  ! Proper axis ratio z/x  !
  printmat(41) = sqrt(1.0d0-(def_quantities_real_(39,impt)/def_quantities_real_(37,impt))**2)   ! sqrt(1-(Rz/Rx)^2)   !
  printmat(42) = def_quantities_real_(40,impt)  ! rho_c !
  printmat(43) = def_quantities_real_(41,impt)  ! pre_c  !
  printmat(44) = def_quantities_real_(42,impt)  ! q_c  !
  printmat(45) = def_quantities_real_(43,impt)  ! rho_max !
  printmat(46) = def_quantities_real_(44,impt)  ! pre_max  !
  printmat(47) = def_quantities_real_(45,impt)  ! epsi_max !
  printmat(48) = def_quantities_real_(46,impt)  ! q_max  !
  printmat(49) = def_quantities_real_(47,impt)  ! coord_radius_x_km !
  printmat(50) = def_quantities_real_(48,impt)  ! coord_radius_y_km !
  printmat(51) = def_quantities_real_(49,impt)  ! coord_radius_z_km !
  printmat(52) = def_quantities_real_(50,impt)  ! proper_radius_x_km !
  printmat(53) = def_quantities_real_(51,impt)  ! proper_radius_y_km !
  printmat(54) = def_quantities_real_(52,impt)  ! proper_radius_z_km !
  printmat(55) = def_quantities_real_(53,impt)  ! rho_c_cgs  !
  printmat(56) = def_quantities_real_(54,impt)  ! pre_c_cgs !
  printmat(57) = def_quantities_real_(55,impt)  ! epsi_c_cgs  !
  printmat(58) = def_quantities_real_(56,impt)  ! q_c_cgs  !
  printmat(59) = def_quantities_real_(57,impt)  ! rho_max_cgs  !
  printmat(60) = def_quantities_real_(58,impt)  ! pre_max_cgs !
  printmat(61) = def_quantities_real_(59,impt)  ! epsi_max_cgs  !
  printmat(62) = def_quantities_real_(60,impt)  ! q_max_cgs  !
  printmat(63) = def_quantities_real_(67,impt)  ! dhdr_x  !
  printmat(64) = def_quantities_real_(68,impt)  ! dhdr_y  !
  printmat(65) = def_quantities_real_(69,impt)  ! dhdr_z  !
  printmat(66) = def_quantities_real_(121,impt) ! chi_cusp  !

  open(100,file=filename,status='unknown')
  do i=1,66
    write(100,'(1p,1e23.15)', ADVANCE="NO")   printmat(i)
  end do
  write(100,*) ""
  close(100)
!
  deallocate(printmat)
!
end subroutine printout_physq_BNS_plot_mpt
