subroutine read_parameter_drot
  use def_matter_parameter
  implicit none
  open(1,file='rnspar_drot.dat',status='old')
  read(1,'(2x,a2,10x,1p,e14.6)') ROT_LAW, DRAT_A2DR
  if (ROT_LAW.eq."DR") then
    read(1,'(1p,2e14.6)') index_DR, A2DR
  end if
  if (ROT_LAW.eq."OJ") then
    read(1,'(1p,2e14.6)') index_DRq, A2DR
    read(1,'(1p,2e14.6)') index_DRp, B2DR
    read(1,'(1p,2e14.6)') DRAT_A2DR, DRAT_B2DR
    if (index_DRq.ne.3.0) stop 'NOT KEPLER'
  end if
  close(1)
end subroutine read_parameter_drot
