subroutine write_parameter_drot
  use def_matter_parameter
  implicit none
  character(len=28) dum
  character(len=50) char_50(4)
!
  open(1,file='rnspar_drot.dat',status='old')
  read(1,'(a28,a50)') dum, char_50(1)  ! For ROT_LAW, DROT_ratio
  if (ROT_LAW.eq."DR") then
    read(1,'(a28,a50)') dum, char_50(2)  ! For index_DR, A2DR
  end if
  if (ROT_LAW.eq."OJ") then
    read(1,'(a28,a50)') dum, char_50(2)  ! For index_DR, A2DR
    read(1,'(a28,a50)') dum, char_50(3)  ! For index_DRp, B2DR
    read(1,'(a28,a50)') dum, char_50(4)  ! For DRAT_A2DR, DRAT_B2DR
  end if
  close(1)
!
  open(1,file='rnspar_drot.las',status='unknown')
  write(1,'(2x,a2,10x,1p,e14.6,a50)') ROT_LAW, DRAT_A2DR, char_50(1)
  if (ROT_LAW.eq."DR") then
    write(1,'(1p,2e14.6,a50)')      index_DR, A2DR,       char_50(2)
  end if
  if (ROT_LAW.eq."OJ") then
    write(1,'(1p,2e14.6,a50)')      index_DRq, A2DR,      char_50(2)
    write(1,'(1p,2e14.6,a50)')      index_DRp, B2DR,      char_50(3)
    write(1,'(1p,2e14.6,a50)')      DRAT_A2DR, DRAT_B2DR, char_50(4)
  end if
  close(1)
end subroutine write_parameter_drot
