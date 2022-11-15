subroutine subio_fluid_1D(istat,iseq,char)
!
  use phys_constant, only : pi
!  use CB_fR_mesh_grav
  use CB_fR_param_iomis
  use grid_parameter_1D	!nrf, nrg
  use def_matter_1D		!emd
  use coordinate_grav_r_1D, only : rg
  implicit none
!
  real(8) :: rftmp(0:nnrg), emddiff, dr0, rmeps
  integer :: ir, nrftmp, nrout, iddd, iwrite
  integer, intent(in) :: istat, iseq
!
  character*4 char
  character*3 moji
!
! ------------------------------------------------------------
!     istat = 0 and iseq 0 -> open ini file and close it 
!                             and open las file
!     istat = 2            -> close las file and open nxt file
! ------------------------------------------------------------
!
  moji = 'nxt'
  if (istat == 0.and.iseq == 0) moji = 'ini'
  if (istat == 1)               moji = 'las'
  if (istat == 2)               moji = 'nxt'
  if (istat == 1) then
    open(2,file='rnsflu_1D.'//moji,status='unknown')
  else 
    open(2,file='rnsflu_1D.'//moji,status='old')
  end if
!
! -------------
! --- Test. ---
! -------------
!
  if (istat == -1) then
!
    emd(0) = emdc
    emd(1:nrf) = emdc*dsin(pi*rg(1:nrf))/(pi*rg(1:nrf))
    ber = 0.8d0
    radi = 0.9d0
!
    return
  end if
!
! -------------
! --- Read. ---
! -------------
!
  if (istat == 0) then
!
! --- Fluid variables
!
    read(2,'(5i5)') nrftmp
    do ir = 0, nrftmp
      read(2,'(6es20.12)') rftmp(ir), emd(ir)
    end do
    read(2,'(6es20.12)') ber, radi
!
    if (nrftmp /= nrf) call interpf_ini(nrftmp,rftmp,emd)
!
  end if
!
! --------------
! --- Write. ---
! --------------
!
  if (istat /= 0) then
!!!
    nrout = nrg - nrgin
    emddiff = 0.0d0
    dr0 = rgmid/dble(nrgin)
    rmeps = 3.0d-05
    iddd = 0
    iwrite = 0
!!!
!
! --- Fluid variables
!
    write(2,'(5i5)') nrf
    do ir = 0, nrf
      write(2,'(6es20.12)') rg(ir), emd(ir)
    end do
    write(2,'(6es20.12)') ber, radi
!
    write(2,'(5i5)') nrg, nrf 
    write(2,'(1i5, 2es10.3)') nrout, dr0, rgout
    write(2,'(1i5, 2es10.3)') iter_max, eps, conv_ini
    write(2,'(2es10.3)') conv0_gra,conv0_den
    write(2,'(3i5, 3x, a4)') itype, iwrite, iddd, char
!
    write(2,'(5i5)') numseq
    write(2,'(2es14.6)') restmass_sph, gravmass_sph
    write(2,'(2es14.6)') emdc, pinx
    write(2,'(2es14.6)') emddiff, rmeps
!
  end if
!
  close(2)
!
end subroutine subio_fluid_1D
