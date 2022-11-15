subroutine qeos_initialize
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_qeos_parameter   !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
!
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  integer :: ii, iphase
!
  abc = 0.0d0
  abi = 0.0d0
  abccgs = 0.0d0
  abcene = 0.0d0
  abiene = 0.0d0
  abch = 0.0d0
  abih = 0.0d0
  rhosurf_cgs = 2.0d0*2.67d14     ! g/cm^(-3)

  open(850,file='qeos_parameter.dat',status='old')
  read(850,'(2es22.15)')    rhosurf_cgs, eneconst_cgs
  read(850,'(17x,1i5,es22.15)') nphase, rhoini_cgs
  do ii = 1, nphase 
    read(850,'(2es22.15)') abccgs(ii), abi(ii)
  end do
!
  close(850)
!
! --  cgs to g = c = msol = 1 unit.
! --  assume pre = pre_0 dyn/cm^2 at rho = rho_0 gr/cm^3.
! --  typically pre_0 = 1.0d+37 dyn/cm^2 
! --  and       rho_0 = 1.0d+16  gr/cm^3.
! --  rescale interface values
!
  facrho = (g/c**2)**3*solmas**2
  facpre = g**3*solmas**2/c**8
  rhosurf_cgs = rhosurf_cgs*2.67d14  
!
!
!  
  do ii=1, nphase  
    abc(ii) = facpre/facrho**abi(ii)*abccgs(ii)
  end do
!
  do ii=1, nphase
    abcene(ii) = abc(ii)/(abi(ii)-1.0d0)
    abiene(ii) = abi(ii)
    abch(ii) = abc(ii)
    abih(ii) = abi(ii)-1.0d0
    abchdot(ii)= abih(ii)*abch(ii)
    abihdot(ii)= abih(ii) - 1.0d0

  end do
  do ii=nphase+1, 2*nphase
    abch(ii) = abcene(ii-nphase)
    abih(ii) = abiene(ii-nphase)-1.0d0
    abchdot(ii)= abih(ii)*abch(ii)
    abihdot(ii)= abih(ii) - 1.0d0
  end do

! calculate abcene abiene abch abih
  open(860,file='peos_parameter_output.dat',status='unknown')
  write(860,'(a1,8x,i5)')'#', nphase
  do ii = 1, nphase
    write(860,'(i5,10es13.5)') ii, abc(ii), abi(ii), abccgs(ii) 
  end do
  close(860)
!
  rhoini_gcm1 = facrho*rhoini_cgs
  rhosurf_gcm1 = facrho*rhosurf_cgs
  eneconst_gcm1 = eneconst_cgs/c**2
!  emdini_gcm1 = abc(iphase)*rhoini_gcm1**(abi(iphase)-1.0d0)
!  do ii=1, 2*nphase
!    write(6,*) 'hdot                               h'
!    write(6,*) abchdot(ii), abihdot(ii), abch(ii), abih(ii)
!  end do
! 
end subroutine qeos_initialize
