subroutine peos_initialize
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter   !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
!
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  real(8) :: eneqc, preqc
  integer :: ii, iphase
!
  open(850,file='peos_parameter.dat',status='old')
!
  read(850,'(8x,1i5,es13.5)') nphase, rhoini_cgs
  read(850,'(2es13.5)') rho_0, pre_0
  read(850,'(2es13.5)') rhocgs(nphase+1), sgma
  do ii = nphase, 0, -1
    read(850,'(2es13.5)') rhocgs(ii), abi(ii)
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
!
  do ii = 0, nphase
    rhoi(ii) = facrho*rhocgs(ii)
  end do
!
  call peos_lookup(rho_0,rhocgs,iphase)
!
  write(6,*) "Reference density rho0 at iphase = ", iphase
!    
  abc(iphase) = pre_0/rho_0**abi(iphase)
  abc(iphase) = facpre/facrho**abi(iphase)*abc(iphase)
  abccgs(iphase) = pre_0/(rho_0**abi(iphase))
!
  if (iphase.gt.0) then
    do ii = iphase-1, 0, -1
      abc(   ii) = rhoi(  ii)**(abi(ii+1)-abi(ii))*abc(   ii+1)
      abccgs(ii) = rhocgs(ii)**(abi(ii+1)-abi(ii))*abccgs(ii+1)
    end do
  end if
  if (iphase.lt.nphase) then
    do ii = iphase+1, nphase
      abc(   ii) = rhoi(  ii-1)**(abi(ii-1)-abi(ii))*abc(   ii-1)
      abccgs(ii) = rhocgs(ii-1)**(abi(ii-1)-abi(ii))*abccgs(ii-1)
    end do
  end if
!
  do ii = 0, nphase
    qi(ii) = abc(ii)*rhoi(ii)**(abi(ii)-1.0d0)
  end do
!
  hi(0) = 1.0d0
  do ii = 1, nphase
    fac2 = abi(ii)/(abi(ii) - 1.0d0)
    hi(ii) = hi(ii-1) + fac2*(qi(ii) - qi(ii-1))
  end do
!
!================================= Quark core =================================
!
  cbar           = hi(nphase)/rhoi(nphase)**sgma  ! continuity of h at rho(nphase)

  eneqc          = (hi(nphase) - qi(nphase))*rhoi(nphase)   ! h = ene/rho + q

  preqc          = qi(nphase)*rhoi(nphase) 

  constqc        = preqc - sgma*eneqc

  rhoi(nphase+1) = facrho*rhocgs(nphase+1)

  hi(nphase+1)   = cbar*rhoi(nphase+1)**sgma

  qi(nphase+1)   = (sgma*cbar*rhoi(nphase+1)**sgma + &
                 &  constqc/rhoi(nphase+1))/(sgma+1.0d0)          
!       
!==============================================================================
!
  open(860,file='peos_parameter_output.dat',status='unknown')
  write(860,'(a1,8x,i5,2es13.5)')'#', nphase, cbar, sgma
  write(860,*) "i, abc, abi, rhoi, qi, hi, prei, enei, abccgs, rhocgs, precgs "
  do ii = 0, nphase
    write(860,'(i5,15es13.5)') ii, abc(ii), abi(ii), rhoi(ii), &
    &    qi(ii), hi(ii), qi(ii)*rhoi(ii), (hi(ii)-qi(ii))*rhoi(ii),  &
    &    abccgs(ii), rhocgs(ii), abccgs(ii)*rhocgs(ii)**abi(ii)
  end do
  ii = nphase+1
  write(860,'(i5,15es13.5)') ii, abc(ii), sgma, rhoi(ii), &
    &    qi(ii), hi(ii), qi(ii)*rhoi(ii), (hi(ii)-qi(ii))*rhoi(ii),  &
    &    abccgs(ii), rhocgs(ii), qi(ii)*rhoi(ii)/facpre
  close(860)
!
  rhoini_gcm1 = facrho*rhoini_cgs
  call peos_lookup(rhoini_gcm1,rhoi,iphase)

  if(iphase==(nphase+1)) then
    emdini_gcm1 = (sgma*cbar*rhoini_gcm1**(sgma+1.0d0) + constqc)/(sgma+1.0d0)/rhoini_gcm1
  else
    emdini_gcm1 = abc(iphase)*rhoini_gcm1**(abi(iphase)-1.0d0)
  end if
  write(6,*) "rhoini iphase = ", iphase
  write(6,*) "rhoini_cgs    = ", rhoini_cgs
  write(6,*) "rhoini_gcm1   = ", rhoini_gcm1
  write(6,*) "emdini_gcm1   = ", emdini_gcm1
!
end subroutine peos_initialize


subroutine peos_initialize_cactus(dir_path)
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter   !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
  character*400, intent(in) :: dir_path
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  real(8) :: eneqc, preqc
  integer :: ii, iphase
!
  open(850,file=trim(dir_path)//'/'//'peos_parameter.dat',status='old')
  read(850,'(8x,1i5,es13.5)') nphase, rhoini_cgs
  read(850,'(2es13.5)') rho_0, pre_0
  read(850,'(2es13.5)') rhocgs(nphase+1), sgma
  do ii = nphase, 0, -1
    read(850,'(2es13.5)') rhocgs(ii), abi(ii)
  end do
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
!
  do ii = 0, nphase
    rhoi(ii) = facrho*rhocgs(ii)
  end do
!
  call peos_lookup(rho_0,rhocgs,iphase)
!    
  abc(iphase) = pre_0/rho_0**abi(iphase)
  abc(iphase) = facpre/facrho**abi(iphase)*abc(iphase)
  abccgs(iphase) = pre_0/(rho_0**abi(iphase))
!
  if (iphase.gt.0) then
    do ii = iphase-1, 0, -1
      abc(   ii) = rhoi(  ii)**(abi(ii+1)-abi(ii))*abc(   ii+1)
      abccgs(ii) = rhocgs(ii)**(abi(ii+1)-abi(ii))*abccgs(ii+1)
    end do
  end if
  if (iphase.lt.nphase) then
    do ii = iphase+1, nphase
      abc(   ii) = rhoi(  ii-1)**(abi(ii-1)-abi(ii))*abc(   ii-1)
      abccgs(ii) = rhocgs(ii-1)**(abi(ii-1)-abi(ii))*abccgs(ii-1)
    end do
  end if
!
  do ii = 0, nphase
    qi(ii) = abc(ii)*rhoi(ii)**(abi(ii)-1.0d0)
  end do
!
  hi(0) = 1.0d0
  do ii = 1, nphase
    fac2 = abi(ii)/(abi(ii) - 1.0d0)
    hi(ii) = hi(ii-1) + fac2*(qi(ii) - qi(ii-1))
  end do
!
!================================= Quark core =================================
!
  cbar           = hi(nphase)/rhoi(nphase)**sgma  ! continuity of h at rho(nphase)

  eneqc          = (hi(nphase) - qi(nphase))*rhoi(nphase)   ! h = ene/rho + q

  preqc          = qi(nphase)*rhoi(nphase)

  constqc        = preqc - sgma*eneqc

  rhoi(nphase+1) = facrho*rhocgs(nphase+1)

  hi(nphase+1)   = cbar*rhoi(nphase+1)**sgma

  qi(nphase+1)   = (sgma*cbar*rhoi(nphase+1)**sgma + &
                 &  constqc/rhoi(nphase+1))/(sgma+1.0d0)
!       
!==============================================================================
!
!
!  open(860,file='peos_parameter_output.dat',status='unknown')
!  write(860,'(a1,8x,i5,2es13.5)')'#', nphase, cbar, sgma
!  do ii = 0, nphase
!    write(860,'(i5,15es13.5)') ii, abc(ii), abi(ii), rhoi(ii), &
!    &    qi(ii), hi(ii), qi(ii)*rhoi(ii), (hi(ii)-qi(ii))*rhoi(ii),  & 
!    &    abccgs(ii), rhocgs(ii), abccgs(ii)*rhocgs(ii)**abi(ii)
!  end do
!  ii = nphase+1
!  write(860,'(i5,15es13.5)') ii, abc(ii), sgma, rhoi(ii), &
!    &    qi(ii), hi(ii), qi(ii)*rhoi(ii), (hi(ii)-qi(ii))*rhoi(ii),  &
!    &    abccgs(ii), rhocgs(ii), qi(ii)*rhoi(ii)/facpre
!  close(860)
!
  rhoini_gcm1 = facrho*rhoini_cgs
  call peos_lookup(rhoini_gcm1,rhoi,iphase)

  if(iphase==(nphase+1)) then
    emdini_gcm1 = (sgma*cbar*rhoini_gcm1**(sgma+1.0d0) + constqc)/(sgma+1.0d0)/rhoini_gcm1
  else
    emdini_gcm1 = abc(iphase)*rhoini_gcm1**(abi(iphase)-1.0d0)
  end if

!
end subroutine peos_initialize_cactus
