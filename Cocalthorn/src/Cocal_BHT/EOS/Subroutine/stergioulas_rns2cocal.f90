subroutine stergioulas_rns2cocal
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_teos_parameter   !rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
!
  integer :: nlines
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  integer :: ii, iphase
  real(8) :: hh, pre, rho, ene, mb
!
  facrho = (g/c**2)**3*solmas**2
  facpre = g**3*solmas**2/c**8
  mb = 1.66d-24
!
  open(850,file='EOS.dat',status='unknown')
  read(850,*) nlines
  do ii=0,nlines-1
!  e/c^2 [g/cm^3],  p [dyn/cm^2],  c^2 log_e(h/c^2) [cm^2/s^2],  n_B [cm^{-3}].
    read(850,*) enei(ii), prei(ii), hi(ii), rhoi(ii)
  end do
  close(850)
!
  open(860,file='teos_parameter.dat',status='unknown')
  do ii=0,nlines-1
    rhoi(ii) = rhoi(ii)*mb
    hh = hi(ii)/c**2
    hi(ii) = dexp(hh)
    write(860,'(1p,7e23.15)') enei(ii), prei(ii), hi(ii), rhoi(ii)
  end do
  close(860)
!
end subroutine stergioulas_rns2cocal
