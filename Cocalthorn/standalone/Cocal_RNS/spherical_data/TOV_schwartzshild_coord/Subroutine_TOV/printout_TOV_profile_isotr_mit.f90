subroutine printout_TOV_profile_isotr_mit(file_ctl)
  use def_quantities
!  use def_qeos_parameter, only : enesurf_gcm1, aq
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  implicit none
  character(2) :: file_ctl
!
  if(file_ctl.eq.'op')  open (8, file='ovlas.dat', status='unknown')
!
  ene  = yn_iso(2)
  call mit_ene2prho(ene,pre,rho0)
  rhocgs=0.0d0
  precgs = pre*c**8/(g**3*solmas**2)
  enecgs = ene*c**8/(g**3*solmas**2)
  xecgs = xe*g*solmas*1.0d-5/c**2
  write(8, '(17(es15.7))') xecgs, yn_iso(1), yn_iso(2), yn_iso(3), yn_iso(4), &
  &     yn_iso(5), yn_iso(6), yn_iso(8), yn_iso(9), rho0, rhocgs, pre, &
  &                        precgs, ene, enecgs, xe, xe/radi
! ene= rho*(1+epsilon)
  if(file_ctl.eq.'cl')  close(8)
!
end subroutine printout_TOV_profile_isotr_mit
