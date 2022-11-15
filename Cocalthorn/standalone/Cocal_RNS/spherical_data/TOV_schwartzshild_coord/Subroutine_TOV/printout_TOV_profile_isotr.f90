subroutine printout_TOV_profile_isotr(file_ctl)
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  implicit none
  character(2) :: file_ctl
  real(8) :: qq
!
  if(file_ctl.eq.'op')  open (8, file='ovlas.dat', status='unknown')
!  open(8, file="ovlas.dat", status="old", position="append")
!
  call peos_h2qprho(yn_iso(2), qq, pre, rho0, ene)

!  if(file_ctl.eq.'cl')   then
!    write(6,*)"*********************************************************************"
!    write(6,'(1p,6e23.15)')  yn_iso(2), qq, pre, rho0, ene
!  end if

  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  xecgs = xe*g*solmas*1.0d-5/c**2
  write(8, '(17(es15.7))') xecgs, yn_iso(1), yn_iso(2), yn_iso(3), yn_iso(4), &
  &     yn_iso(5), yn_iso(6), yn_iso(8), yn_iso(9), rho0, rhocgs, pre, &
  &                        precgs, ene, epsiloncgs, xe, xe/radi
!
!  close(8)
  if(file_ctl.eq.'cl')  close(8)
!
end subroutine printout_TOV_profile_isotr
