subroutine read_parameter_bh_mpt(impt)
  use phys_constant, only : long, pi
  use grid_parameter, only : num_sol_seq
  use def_bh_parameter
  use def_binary_parameter, only : mass_ratio
  implicit none
  integer,intent(in)  :: impt
  character(len=1) :: np(5) = (/'1', '2','3', '4', '5'/)
  open(1,file='bbhpar_mpt'//np(impt)//'.dat',status='old')
  read(1,'(1p,2e14.6)') ome_bh, spin_bh
  read(1,'(1p,2e14.6)') alph_bh, psi_bh
  read(1,'(1p,2e14.6)') th_spin_bh_deg, phi_spin_bh_deg
  read(1,'(2(2x,a2),6x,1p,e14.6)') bh_bctype, bh_soltype, mass_ratio
  read(1,'(1p,e14.6,i5)') rgin_deform, num_sol_seq
  close(1)
  th_spin_bh  = pi*th_spin_bh_deg/180.0d0
  phi_spin_bh = pi*phi_spin_bh_deg/180.0d0
  bh_sptype = 'SP'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq.  0.0) bh_sptype = 'Xp'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq. 90.0) bh_sptype = 'Yp'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq.180.0) bh_sptype = 'Xm'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq.270.0) bh_sptype = 'Ym'
  if (th_spin_bh_deg.eq.  0.0) bh_sptype = 'Zp'
  if (th_spin_bh_deg.eq.180.0) bh_sptype = 'Zm'
end subroutine read_parameter_bh_mpt
