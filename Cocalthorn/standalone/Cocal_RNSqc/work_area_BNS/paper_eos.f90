PROGRAM paper_eos
  IMPLICIT NONE
!
  real(8) :: rho_0, pre_0, facrho, facpre, fac2, gg, cc, ss
  integer :: j,ii, iphase, nphase
  real(8) :: rhoini_cgs, rhocgs(0:100), abi(0:100)

  character(LEN=100) :: inputfile, cmd, tempchar1, tempchar2
  character(len=5) :: chy,chx
!
  open(850,file='peos_parameter.dat',status='old')
  read(850,'(8x,1i5,es13.5)') nphase, rhoini_cgs
  read(850,'(2es13.5)') rho_0, pre_0
  do ii = nphase, 0, -1
    read(850,'(2es13.5)') rhocgs(ii), abi(ii)
  end do
  close(850)

  write(6,*) "$\ \log(P_2)\ $ & ", "$\ \Gamma_1\ $ & ", "$\ \Gamma_2\ $ & ", "$\ \Gamma_3\ $ & ", &
   &   "$\ \Gamma_4\ $ & ", "$\ \log((\GR_0)_1)\ $ & ", "$\ \log((\GR_0)_2)\ $ & ", "$\ \log((\GR_0)_3)\ $ & "

  write(6,'(f7.3,a3)', ADVANCE = "NO") log10(pre_0), " & "
  do j=1,nphase
    write(6,'(f6.3,a3)', ADVANCE = "NO") abi(j), " & " 
  end do
  do j=1,nphase-1
    write(6,'(f7.3,a3)', ADVANCE = "NO") log10(rhocgs(j)), " & "
  end do
  write(6,*) ""

END PROGRAM paper_eos

