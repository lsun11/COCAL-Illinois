subroutine peos2teos(nlines)
!
  use phys_constant        !g,c,solmas,nnpeos
  use def_peos_parameter   !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  implicit none
!
  integer, intent(in) :: nlines
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss, mb
  integer :: ii, iphase
  real(8) :: dlogrho, rhoteos, preteos, qteos, hh, pre, rho, ene
  real(8) :: ene1,pre1,h1,rho1,q1, ide, hprev
!
  mb = 1.66d-24

  open(850,file='peos_parameter.dat',status='old')
!
  read(850,'(8x,1i5,es13.5)') nphase, rhoini_cgs
  read(850,'(2es13.5)') rho_0, pre_0
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
  open(860,file='peos_parameter_output.dat',status='unknown')
  write(860,'(a1,8x,i5)')'#', nphase
  do ii = 0, nphase
    write(860,'(i5,10es13.5)') ii, abc(ii), abi(ii), rhoi(ii), &
    &                           qi(ii), hi(ii), abccgs(ii), rhocgs(ii), &
    &                           abccgs(ii)*rhocgs(ii)**abi(ii)
  end do
  close(860)
!
  rhoini_gcm1 = facrho*rhoini_cgs
  call peos_lookup(rhoini_gcm1,rhoi,iphase)
  emdini_gcm1 = abc(iphase)*rhoini_gcm1**(abi(iphase)-1.0d0)
!
!
! Check for a certain baryon number density 
!
!  rhoteos = mb*3.469614697607567d+40
!  call peos_lookup(rhoteos,rhocgs,iphase)
!  preteos = abccgs(iphase)*rhoteos**abi(iphase)
!  qteos   = preteos/rhoteos/(c**2)
!  call peos_q2hprho(qteos, hh, pre, rho, ene)
!  ene = ene/facrho
!  rho = rho/facrho
!  pre = pre/facpre
!  write(6,'(1p,7e23.15)')  ene, pre, dlog(hh)*c**2, rho/mb
!stop
!
  open(870,file='peos_parameter_ephrho.dat',status='unknown')
  write(870,'(a1,8x,i5)')'#', nphase
  do ii = 0, nphase
    pre1 = qi(ii)*rhoi(ii)
    write(870,'(1p,10e23.15)') rhoi(ii)*hi(ii)-pre1, pre1, hi(ii), rhoi(ii)
  end do
  write(870,'(a1)')'#'
  do ii = 0, nphase
    pre1 = qi(ii)*rhoi(ii)
    write(870,'(1p,10e23.15)') dlog10(rhoi(ii)*hi(ii)-pre1), dlog10(pre1), &
       &  dlog10(hi(ii)), dlog10(rhoi(ii))
  end do
  close(870)

  hprev=-1.0d0
  write(6,*) "rhocgs(nphase)=", rhocgs(nphase), "     nlines=", nlines
  ide = -0.2
  dlogrho = (dlog10(rhocgs(nphase))-ide)/dble(nlines)
  dlogrho = 18.0/dble(nlines)
!  write(6,*)  "dlogrho=", dlogrho
  open(850,file='teos_parameter.dat',status='unknown')
  open(860,file='teos_parameter_cgs.dat',status='unknown')
!  do ii=1,nlines
  do ii=1,6000
!    rhoteos = (10.0d0**ide)*(10.0d0**(dble(ii)*dlogrho))
    if(ii<=500)  then
      rhoteos = dble(ii)*0.2d0*10.0d0**10
    else if (ii>500 .and. ii<=1000)  then
      rhoteos = (1.0d0 + dble(ii-500)*0.2d0)*10.0d0**12
    else
      rhoteos = (1.0d0 + dble(ii-1000)*0.2d0)*10.0d0**14
    end if
    call peos_lookup(rhoteos,rhocgs,iphase)
    preteos = abccgs(iphase)*rhoteos**abi(iphase)

    rho1 = rhoteos*facrho
    pre1 = preteos*facpre
    q1   = pre1/rho1

    call peos_q2hprho(q1, h1, pre1, rho1, ene1)

    ene = ene1/facrho 
    rho = rho1/facrho   
    pre = pre1/facpre   
    hh = (ene*c**2+pre)/rho

    if (h1 .ne. hprev)  then
      write(850,'(1p,4e23.15)')  ene1, pre1, h1, rho1
      write(860,'(1p,4e23.15)')  ene,  pre,  hh, rho
    end if

    hprev=h1
!    rhoi(ii) = rhocgs(ii)*facrho
!    prei(ii) = precgs(ii)*facpre
!    qi(ii)   = prei(ii)/rhoi(ii)
  end do
  close(850)
  close(860)

end subroutine peos2teos
