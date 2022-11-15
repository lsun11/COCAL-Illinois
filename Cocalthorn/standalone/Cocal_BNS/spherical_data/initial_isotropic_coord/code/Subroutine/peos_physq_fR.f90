subroutine peos_physq_fR(istat,iseq)
!
  use CB_fR_param_flphy, only : ahoadm, ahores
  use def_matter_1D, only : radi, ber
  use grid_parameter_1D, only : nrf, chrot, chgra, chope, nrg, emdc, &
  &                          gravmass_sph, pinx, restmass_sph, &
  &                          rgout
  use weight_grav_1D, only : wgdrg
  implicit none
  real(8) :: admmas, ahokom, fitadm, fitkom, fitvir, pregc, rhogc, &
  &          rkomas, virial, hhc, ene
  integer, intent(in) :: istat, iseq
!
  if (iseq == 1) then
  open(10,file='frgphy.dat',status='old')
  open(11,file='frgphydp.dat',status='old')
  open(31,file='frgphyplot.dat',status='old')
  open(34,file='frgasym.dat',status='old')
  end if
  if (istat == 2) write(10,*) ' == Solution does not converge. == '
!
!testtest
  ahokom = 0.0d0
  virial = 0.0d0
  admmas = 0.0d0
  rkomas = 0.0d0
  fitvir = 0.0d0
  fitadm = 0.0d0
  fitkom = 0.0d0

!testtest
!
  write(10,'(2i7, 6es14.6)') nrf, nrg, pinx, rgout
  write(10,'(6es14.6)') radi, ber, emdc
  write(10,'(6es14.6)') fitvir, fitadm, fitkom
  write(10,'(6es14.6)') ahores, restmass_sph
  write(10,'(6es14.6)') ahoadm, gravmass_sph
!
  write(6,'(2i7, 6es14.6)') nrf, nrg, pinx, rgout
  write(6,'(6es14.6)') radi, ber, emdc
  write(6,'(6es14.6)') fitvir, fitadm, fitkom
  write(6,'(6es14.6)') ahores, restmass_sph
  write(6,'(6es14.6)') ahoadm, gravmass_sph
!
  write(11,'(2i7, 6es14.6)') nrf, nrg, pinx
  write(11,'(3es23.15)') radi, ber
  write(11,'(3es23.15)') emdc
  write(11,'(3es23.15)') ahores, ahoadm, ahokom 
  write(11,'(3es23.15)') virial, admmas, rkomas
  write(11,'(3es23.15)') fitvir, fitadm, fitkom
  write(11,'(3es23.15)') restmass_sph, gravmass_sph
!
!  rhogc = emdc**pinx
!  pregc = emdc*rhogc
  call peos_q2hprho(emdc, hhc, pregc, rhogc, ene)
  write(31,'(40es16.8)') radi, emdc, rhogc, pregc, &
     &               ahores, ahoadm, ahokom, &
     &               virial, admmas, rkomas, &
     &               fitvir, fitadm, fitkom, restmass_sph, gravmass_sph
!
end subroutine peos_physq_fR
