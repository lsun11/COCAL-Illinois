subroutine peos_initialize_mpt(impt)
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter   !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
  integer,intent(in)  :: impt
  character(len=1) :: np(2) = (/'1', '2'/)
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  real(8) :: ene, pre, hh, rho
  integer :: ii, iphase, iofile, idum, i0
  real(8), external      :: lagint_4th, lagint_2nd
  real(8)                :: x4(4), f4(4)
!
! --  cgs to g = c = msol = 1 unit.
! --  assume pre = pre_0 dyn/cm^2 at rho = rho_0 gr/cm^3.
! --  typically pre_0 = 1.0d+37 dyn/cm^2
! --  and       rho_0 = 1.0d+16  gr/cm^3.
!
  facrho = (g/c**2)**3*solmas**2
  facpre = g**3*solmas**2/c**8

!  open(850,file='teos_parameter.dat',status='old')
!  nphase=0
!  do
!   read(850,*,iostat=iofile)  rho_0
!   if (iofile > 0)  then
!     write(6,*) "Problem reading file teos_parameter.dat...exiting"
!     stop
!   else if (iofile < 0) THEN
!!    nphase is the total number of phases. The total number of points is nphase+1
!     nphase = nphase - 1
!     write(6,*) "Total number of lines is  =", nphase+1
!     write(6,*) "Total number of phases is =", nphase
!     exit
!   else
!     nphase = nphase + 1
!   end if
!!   write(6,*) nphase, rho_0
!  end do
!  close(850)

  open(850,file='teos_parameter_mpt'//np(impt)//'.dat',status='old')
! nphase is the total number of phases. The total number of points is nphase+1
  read(850,'(i5)') idum   ! number of points (lines)
  nphase=idum-1
  do ii=0,nphase
    read(850,'(1p,4e23.15)') ene, pre, hh, rho
    if (ii==0 .or. ii==nphase)   write(6,'(i5,1p,4e23.15)') ii, ene, pre, hh, rho
    enei(ii) = ene
    rhoi(ii) = rho
    prei(ii) = pre
    hi(ii)   = hh
    qi(ii)   = prei(ii)/rhoi(ii) 
!    rhoi(ii) = rhocgs(ii)*facrho
!    prei(ii) = precgs(ii)*facpre
!    write(6,'(i4,4e20.12)') ii, rhoi(ii), prei(ii), hi(ii), qi(ii)
  end do
  close(850)
  write(6,*) "i, rhoi, prei, hi, qi"
  write(6,'(i5,1p,4e23.15)') 0, rhoi(0), prei(0), hi(0), qi(0)
  write(6,'(i5,1p,4e23.15)') nphase, rhoi(nphase), prei(nphase), hi(nphase), qi(nphase)
!
  open(850,file='peos_parameter_mpt'//np(impt)//'.dat',status='old')
  read(850,'(8x,1i5,es13.5)') idum, rhoini_cgs
  read(850,'(2es13.5)') kappa_crust, gamma_crust
  close(850)
  write(6,'(a21,1p,e23.15)') "******* rhoini_cgs = ", rhoini_cgs
  rhoini_gcm1 = facrho*rhoini_cgs

!  rhoini_gcm1 = 1.73d-03
!  rhoini_cgs  = rhoini_gcm1/facrho

  call peos_lookup(rhoini_gcm1,rhoi,iphase)
  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = rhoi(i0:i0+3)
  f4(1:4) = qi(i0:i0+3)
  emdini_gcm1 = lagint_4th(x4,f4,rhoini_gcm1)

!  emdini_gcm1 = 1.132621660533359d-01  !FPS
!  emdini_gcm1 = 1.600980540676673d-01  !gamma=2
  write(6,'(a26,1p,2e23.15)') "(rhoini_gcm1,emdini_gcm1)=", rhoini_gcm1, emdini_gcm1

end subroutine peos_initialize_mpt

subroutine peos_initialize_mpt_cactus(impt, dir_path)
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter   !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
  character*400, intent(in) :: dir_path
  integer,intent(in)  :: impt
  character(len=1) :: np(2) = (/'1', '2'/)
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  real(8) :: ene, pre, hh, rho
  integer :: ii, iphase, iofile, idum, i0
  real(8), external      :: lagint_4th, lagint_2nd
  real(8)                :: x4(4), f4(4)
!
! --  cgs to g = c = msol = 1 unit.
! --  assume pre = pre_0 dyn/cm^2 at rho = rho_0 gr/cm^3.
! --  typically pre_0 = 1.0d+37 dyn/cm^2
! --  and       rho_0 = 1.0d+16  gr/cm^3.
!
  facrho = (g/c**2)**3*solmas**2
  facpre = g**3*solmas**2/c**8

!  open(850,file='teos_parameter.dat',status='old')
!  nphase=0
!  do
!   read(850,*,iostat=iofile)  rho_0
!   if (iofile > 0)  then
!     write(6,*) "Problem reading file teos_parameter.dat...exiting"
!     stop
!   else if (iofile < 0) THEN
!!    nphase is the total number of phases. The total number of points is nphase+1
!     nphase = nphase - 1
!     write(6,*) "Total number of lines is  =", nphase+1
!     write(6,*) "Total number of phases is =", nphase
!     exit
!   else
!     nphase = nphase + 1
!   end if
!!   write(6,*) nphase, rho_0
!  end do
!  close(850)

  open(850,file=trim(dir_path)//'/'//'teos_parameter_mpt'//np(impt)//'.dat',status='old')
! nphase is the total number of phases. The total number of points is nphase+1
  read(850,'(i5)') idum   ! number of points (lines)
  nphase=idum-1
  do ii=0,nphase
    read(850,'(1p,4e23.15)') ene, pre, hh, rho
    if (ii==0 .or. ii==nphase)   write(6,'(i5,1p,4e23.15)') ii, ene, pre, hh, rho
    enei(ii) = ene
    rhoi(ii) = rho
    prei(ii) = pre
    hi(ii)   = hh
    qi(ii)   = prei(ii)/rhoi(ii) 
!    rhoi(ii) = rhocgs(ii)*facrho
!    prei(ii) = precgs(ii)*facpre
!    write(6,'(i4,4e20.12)') ii, rhoi(ii), prei(ii), hi(ii), qi(ii)
  end do
  close(850)
  write(6,*) "i, rhoi, prei, hi, qi"
  write(6,'(i5,1p,4e23.15)') 0, rhoi(0), prei(0), hi(0), qi(0)
  write(6,'(i5,1p,4e23.15)') nphase, rhoi(nphase), prei(nphase), hi(nphase), qi(nphase)
!
  open(850,file=trim(dir_path)//'/'//'peos_parameter_mpt'//np(impt)//'.dat',status='old')
  read(850,'(8x,1i5,es13.5)') idum, rhoini_cgs
  read(850,'(2es13.5)') kappa_crust, gamma_crust
  close(850)
  write(6,'(a21,1p,e23.15)') "******* rhoini_cgs = ", rhoini_cgs
  rhoini_gcm1 = facrho*rhoini_cgs

!  rhoini_gcm1 = 1.73d-03
!  rhoini_cgs  = rhoini_gcm1/facrho

  call peos_lookup(rhoini_gcm1,rhoi,iphase)
  i0 = min0(max0(iphase-2,0),nphase-3)
  x4(1:4) = rhoi(i0:i0+3)
  f4(1:4) = qi(i0:i0+3)
  emdini_gcm1 = lagint_4th(x4,f4,rhoini_gcm1)

!  emdini_gcm1 = 1.132621660533359d-01  !FPS
!  emdini_gcm1 = 1.600980540676673d-01  !gamma=2
  write(6,'(a26,1p,2e23.15)') "(rhoini_gcm1,emdini_gcm1)=", rhoini_gcm1, emdini_gcm1

end subroutine peos_initialize_mpt_cactus
