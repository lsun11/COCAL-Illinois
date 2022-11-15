subroutine printout_TOV_profile(file_ctl)
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas
  implicit none
  character(2) :: file_ctl
!
  if(file_ctl.eq.'op')  open (8, file='ovlas.dat', status='unknown')
!
  call peos_h2qprho(yn(2), hh, pre, rho0, ene)
  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  xecgs = xe*g*solmas*1.0d-5/c**2
  write(8, '(13(es15.7))') xecgs, yn(1), yn(2), yn(3), yn(4), &
  &                        yn(5), yn(6), rho0, rhocgs, pre, &
  &                        precgs, ene, epsiloncgs
!
  if(file_ctl.eq.'cl')  close(8)
!
end subroutine printout_TOV_profile
