subroutine read_bht_parameter
  use phys_constant, only : long, pi
  use grid_parameter, only : r_surf, num_sol_seq
  use def_bh_parameter
  use def_bht_parameter
  implicit none

  open(1,file='bhtpar.dat',status='old')
!  read(1,'(i5,4x,a1,3i5)') nrg_1, sw_L1_iter, sw_sepa, sw_quant, sw_spin
  read(1,'(1p,2e23.15)') r_surf, mass_bh
  read(1,'(1p,2e23.15)') ome_bh, spin_bh
  read(1,'(1p,2e23.15)') alph_bh, psi_bh
  read(1,'(1p,2e23.15)') th_spin_bh_deg, phi_spin_bh_deg
  read(1,'(2(2x,a2),15x,1p,e23.15)') bh_bctype, bh_soltype, emdc_ratio
  read(1,'(1p,2e23.15)') l_r_surf, o_r_surf
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
!
end subroutine read_bht_parameter


subroutine read_bht_parameter_cactus(dir_path)
  use phys_constant, only : long, pi
  use grid_parameter, only : r_surf, num_sol_seq
  use def_bh_parameter
  use def_bht_parameter
  implicit none
  character*400, intent(in) :: dir_path

  write(6,*) "BHT directory path:", dir_path

  open(1,file=trim(dir_path)//'/'//'bhtpar.dat',status='old')
!  read(1,'(i5,4x,a1,3i5)') nrg_1, sw_L1_iter, sw_sepa, sw_quant, sw_spin
  read(1,'(1p,2e23.15)') r_surf, mass_bh
  read(1,'(1p,2e23.15)') ome_bh, spin_bh
  read(1,'(1p,2e23.15)') alph_bh, psi_bh
  read(1,'(1p,2e23.15)') th_spin_bh_deg, phi_spin_bh_deg
  read(1,'(2(2x,a2),15x,1p,e23.15)') bh_bctype, bh_soltype, emdc_ratio
  read(1,'(1p,2e23.15)') l_r_surf, o_r_surf
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
!
end subroutine read_bht_parameter_cactus

