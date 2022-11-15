subroutine mit_initialize
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_qeos_parameter   !rhoini_cgs, rhosurf_cgs, aq
  implicit none
!
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  integer :: ii, iphase
!

  open(850,file='mit_parameter.dat',status='old')
  read(850,'(es22.15)') aq
  read(850,'(2es22.15)') rhoini_cgs, rhosurf_cgs
  close(850)

 
!
!
! --  cgs to g = c = msol = 1 unit.
! --  assume pre = pre_0 dyn/cm^2 at rho = rho_0 gr/cm^3.
! --  typically pre_0 = 1.0d+37 dyn/cm^2 
! --  and       rho_0 = 1.0d+16  gr/cm^3.
! --  rescale interface values
!
  facrho = (g/c**2)**3*solmas**2
  facpre = g**3*solmas**2/c**8
  eneini_cgs = rhoini_cgs*c**2
  enesurf_cgs = rhosurf_cgs*c**2
  eneini_gcm1 = eneini_cgs*facpre
  enesurf_gcm1 = enesurf_cgs*facpre
  rhoini_cgs = 0
  rhosurf_cgs = 0  

! calculate abcene abiene abch abih
!
  rhoini_gcm1 = 0
  rhosurf_gcm1 = 0
!  emdini_gcm1 = abc(iphase)*rhoini_gcm1**(abi(iphase)-1.0d0)

 
! 
end subroutine mit_initialize
