subroutine read_parameter_pbh
  use phys_constant, only : long, pi
  use def_bh_parameter
  implicit none
  character(len=3) :: spin_input
  real(long) :: spin_bh_xy
  open(1,file='pbhpar.dat',status='old')
  read(1,'(1p,1e14.6,2x,a3)') mass_pBH, spin_input
  read(1,'(1p,2e14.6)') ome_bh, spin_bh
  read(1,'(1p,2e14.6)') th_spin_bh_deg, phi_spin_bh_deg
  read(1,'(1p,2e14.6)') mom_pBH(1), spin_pBH(1)
  read(1,'(1p,2e14.6)') mom_pBH(2), spin_pBH(2)
  read(1,'(1p,2e14.6)') mom_pBH(3), spin_pBH(3)
  close(1)
  if (spin_input.eq.'ANG') then
    th_spin_bh  = pi* th_spin_bh_deg/180.0d0
    phi_spin_bh = pi*phi_spin_bh_deg/180.0d0
    spin_pBH(1) = spin_bh*dsin(th_spin_bh)*dcos(phi_spin_bh)
    spin_pBH(2) = spin_bh*dsin(th_spin_bh)*dsin(phi_spin_bh)
    spin_pBH(3) = spin_bh*dcos(th_spin_bh)
  else
    spin_bh = dsqrt(spin_pBH(1)**2 + spin_pBH(2)**2 + spin_pBH(3)**2)
    spin_bh_xy = dsqrt(spin_pBH(1)**2 + spin_pBH(2)**2)
    if (spin_bh.eq.0.0d0) then
      th_spin_bh  = 0.0d0 ; th_spin_bh_deg  = 0.0
      phi_spin_bh = 0.0d0 ; phi_spin_bh_deg = 0.0
    else if (spin_bh_xy.eq.0.0d0) then
      th_spin_bh  = dmod(atan2(spin_bh_xy,spin_pbh(3))+2.0d0*pi,2.0d0*pi)
      th_spin_bh_deg  = th_spin_bh*180.0d0/pi
      phi_spin_bh = 0.0d0 ; phi_spin_bh_deg = 0.0
    else
      th_spin_bh  = dmod(atan2(spin_bh_xy,spin_pbh(3))+2.0d0*pi,2.0d0*pi)
      th_spin_bh_deg  = th_spin_bh*180.0d0/pi
      phi_spin_bh = dmod(atan2(spin_pbh(2),spin_pbh(1))+2.0d0*pi,2.0d0*pi)
      phi_spin_bh_deg = phi_spin_bh*180.0d0/pi
    end if
  end if
  bh_sptype = 'SP'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq.  0.0) bh_sptype = 'Xp'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq. 90.0) bh_sptype = 'Yp'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq.180.0) bh_sptype = 'Xm'
  if (th_spin_bh_deg.eq.90.0.or.phi_spin_bh_deg.eq.270.0) bh_sptype = 'Ym'
  if (th_spin_bh_deg.eq.  0.0) bh_sptype = 'Zp'
  if (th_spin_bh_deg.eq.180.0) bh_sptype = 'Zm'
end subroutine read_parameter_pbh
