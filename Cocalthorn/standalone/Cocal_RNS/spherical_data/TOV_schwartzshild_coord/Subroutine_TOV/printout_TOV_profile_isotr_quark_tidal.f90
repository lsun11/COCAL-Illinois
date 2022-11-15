subroutine printout_TOV_profile_isotr_quark_tidal(file_ctl)
  use def_quantities
  use def_TOV_quantities
  use phys_constant, only:g,c,solmas,pi
  implicit none
  character(2) :: file_ctl
  real(8),save :: y_0,ene_0
  real(8) :: c_s, y,co2,co3,co5,deno,rma,rsh,ht
!
  if(file_ctl.eq.'op')  open (8, file='ovlas.dat', status='unknown')
!

  call quark_rho2phenedpdrho(yn_iti(2),pre,hh,ene,dpdrho)
  rho0=yn_iti(2)
!  call peos_h2qprho(yn_iti(2), hh, pre, rho0, ene)
  rhocgs = rho0*(c**6/((g**3)*solmas**2))
  precgs = pre*c**8/(g**3*solmas**2)
  epsiloncgs = ene*c**6/(g**3*solmas**2)
  xecgs = xe*g*solmas*1.0d-5/c**2
  c_s=0
  call quark_sound_speed(c_s,rho0)
  rma = yn_iti(1)
  rsh = yn_iti(9)
  ht = yn_iti(11)

  if (file_ctl.ne.'cl') then
    y  = yn_iti(9)*yn_iti(12)/yn_iti(11)
    y_0 = y
    ene_0 = ene
  else
    y  = y_0 - 4.0d0*pi*ene_0*yn_iti(9)**3/yn_iti(1)
  end if

  compa = yn_iti(1)/yn_iti(9)
  co2=compa*compa
  co3=co2*compa
  co5=co3*co2

  deno = 2.0d0*compa*(6.0d0-3.0d0*y+3.0d0*compa*(5.0d0*y-8.0d0))  &
   &   + 4.0d0*co3*(13.0d0-11.0d0*y+compa*(3.0d0*y-2.0d0)+2.0d0*co2*(1.0d0+y)) &
   &   + 3.0d0*(1.0d0-2.0d0*compa)**2*(2.0d0-y+2.0d0*compa*(y-1.0d0))*dlog(1.0d0-2.0d0*compa)

  k2 = 8.0d0*co5/5.0d0*(1.0d0-2.0d0*compa)**2*(2.0d0+2.0d0*compa*(y-1.0d0)-y)/deno
 
  write(8, '(22(es15.7))') xecgs, yn_iti(1), yn_iti(2), yn_iti(3), yn_iti(4), &
  &     yn_iti(5), yn_iti(6), yn_iti(8), yn_iti(9), hh, rhocgs, pre, &
  &     precgs, ene, epsiloncgs, xe, xe/radi,  yn_iti(11),yn_iti(12),c_s, y, k2
!
  if(file_ctl.eq.'cl')  close(8)
!
end subroutine printout_TOV_profile_isotr_quark_tidal
