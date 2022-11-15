subroutine read_parameter_bh
  use phys_constant, only : long, pi
  use def_bh_parameter
  implicit none
  open(1,file='bbhpar.dat',status='old')
  read(1,'(1p,2e14.6)') ome_bh, spin_bh
  read(1,'(1p,2e14.6)') alph_bh, psi_bh
  read(1,'(1p,2e14.6)') th_spin_bh_deg, phi_spin_bh_deg
  read(1,'(2(2x,a2))') bh_bctype, bh_soltype
  read(1,'(1p,2e14.6)') mass_pBH
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
end subroutine read_parameter_bh
