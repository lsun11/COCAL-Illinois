subroutine IO_output_physq_paper(iseq)
  use def_quantities
  use def_quantities_derived
  implicit none
  integer, intent(in) :: iseq
  real(8) :: pinx_tmp
!
  if (iseq.eq.1) then 
    pinx_tmp = (epsilon_c_normalized/rho_c_normalized - 1.0d0) & 
    &          /p_over_rho_c_normalized
    open(300,file='rnsphy_paper.dat',status='unknown')
    write(300,'("$ n = ",1f3.1,"\quad$   $ M/R = ",1f4.2, &
    &         "\quad$  $ M_0 =",1es12.4,"\quad$  $ M =",1es12.4,"$")') &
    &  pinx_tmp, MoverR_sph, restmass_sph, gravmass_sph
  end if
! 
! Reference for the data 
! Slot # 
! 1  2 ' NS radius along x  = ', coord_radius_x, proper_radius_x
! 3  4 ' Axis ratio y/x     = ', coord_axis_ratio_yx, proper_axis_ratio_yx
! 5  6 ' Axis ratio z/x     = ', coord_axis_ratio_zx, proper_axis_ratio_zx
!    '## Eccentricity in coordinate and Proper Radii ##'
! 7    ' epsilon_max        = ', epsilon_c_normalized
! 8   ' Omega            = ', omega
! 9   ' M_ADM              = ', admmass
! 10   ' Angular momentum J = ', angmom
! 11 ' T/|W|              = ', ToverW
! 12 ' Moment of inertia  = ', I_inertia
! 13 ' surface on z axis  = ', zrb_zp_plus
!
!
write(300,'(3("$",1f6.4," $ $(",1f6.4,")$ & "),   &
&           2("$",1f6.4," $ & "),                 &
&           1("$",1es12.4," $ & "),               &
&           1("$",1es11.3," $ & "),               &
&           1("$",1f6.4," $ & "),                 &
&           1("$",1es11.3," $ & "),               &
&           1("$",1f6.4," $ & "))') &
coord_radius_x, proper_radius_x, &
coord_axis_ratio_yx, proper_axis_ratio_yx, &
coord_axis_ratio_zx, proper_axis_ratio_zx, &
epsilon_c_normalized, omega, &
admmass, angmom, ToverW, I_inertia, zrb_zp_plus
!
end subroutine IO_output_physq_paper
