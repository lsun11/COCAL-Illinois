subroutine subio_1D(istat,iseq,char)
!
  use phys_constant, only : pi
  use def_metric_1D
  use CB_fR_param_iomis
  use grid_parameter_1D	!nrf, nrg
  use def_matter_1D		!ber
  use coordinate_grav_r_1D, only : rg
  implicit none
!
  integer, intent(in) :: istat, iseq
  real(8) :: rgtmp(0:nnrg), dr0, rmeps, emddiff
  integer :: irg, nrgtmp, nrout, iddd, iwrite
!
  character*4 char
  character*3 moji
!
! ------------------------------------------------------------
!     istat = 0 and iseq 0 -> open ini file and close it. 
!     istat = 1            -> open las file and close it.
!     istat = 2            -> open nxt file and close it.
! ------------------------------------------------------------
!
  moji = 'nxt'
  if (istat == 0.and.iseq == 0) moji = 'ini'
  if (istat == 1)               moji = 'las'
  if (istat == 2)               moji = 'nxt'
  if (istat == 1) then
    open(3,file='rnsgra_1D.'//moji,status='unknown')
  else 
    open(3,file='rnsgra_1D.'//moji,status='old')
  end if
!
! -------------
! --- Test. ---
! -------------
!
  if (istat == -1) then
!
    psi(0:nrg) = 1.0d0
    alph(0:nrg) = 1.0d0
    alph(0) = 0.9d0 - emdc
    alph(1:nrf) = 0.9d0 - emdc*dsin(pi*rg(1:nrf))/(pi*rg(1:nrf))
!
  end if
!
! -------------
! --- Read. ---
! -------------
!
  if (istat == 0) then
!
! --- Metric potentials.  
!
  read(3,'(5i5)') nrgtmp
  do irg = 0, nrgtmp
    read(3,'(6es20.12)') rgtmp(irg), psi(irg), alph(irg)
  end do
  read(3,'(6es20.12)') ber, radi
!
  if (nrgtmp /= nrg) then 
  call interpg_ini(nrgtmp,rgtmp,psi)
  call interpg_ini(nrgtmp,rgtmp,alph)
  end if
!
  open(24,file='frgaux.dat',status='old')
  do irg = 0, nrg
    write(24,'(1i5,6es14.6)') irg,rg(irg), psi(irg), alph(irg)
  end do
  close(24)
!        stop
!
  close(3)
!
  end if
!
! --------------
! --- Write. ---
! --------------
!
  if (istat == 1.or.istat == 2) then
!!!
  nrout = nrg - nrgin
  emddiff = 0.0d0
  dr0 = rgmid/dble(nrgin)
  rmeps = 3.0d-05
  iddd = 0
  iwrite = 0
!!!
!
! --- Metric potentials.  
!
  write(3,'(5i5)') nrg
  do irg = 0, nrg
    write(3,'(6es20.12)') rg(irg), psi(irg), alph(irg), emdg(irg)
  end do
  write(3,'(6es20.12)') ber, radi
!
  write(3,'(5i5)') nrg, nrf 
  write(3,'(1i5, 2es10.3)') nrout, dr0, rgout
  write(3,'(1i5, 2es10.3)') iter_max, eps, conv_ini
  write(3,'(2es10.3)') conv0_gra,conv0_den
  write(3,'(3i5, 3x, a4)') itype, iwrite, iddd, char
!
  write(3,'(5i5)') numseq
  write(3,'(2es14.6)') restmass_sph, gravmass_sph
  write(3,'(2es14.6)') emdc, pinx
  write(3,'(2es14.6)') emddiff, rmeps
!
  close(3)
!
! --  For plotter
!
  open(3,file='frggraplot.dat',status='unknown')
  do irg = 0, nrg
    write(3,'(6es20.12)') rg(irg), psi(irg), alph(irg), emdg(irg)
  end do
  close(3)
!
  end if
!
end subroutine subio_1D
