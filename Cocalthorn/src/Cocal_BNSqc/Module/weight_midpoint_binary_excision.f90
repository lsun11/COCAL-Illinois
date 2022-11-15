!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_binary_excision
  use phys_constant, only : long
  implicit none
! weight binary excision
  real(long), pointer :: hwrtpg_ex(:,:,:)
  real(long), pointer :: hwtpg_ex(:,:)
!
contains
!
subroutine allocate_weight_midpoint_binary_excision
  use grid_parameter, only : nrg, ntg, npg 
  use make_array_2d
  use make_array_3d
  implicit none
  call alloc_array3d(hwrtpg_ex,1,nrg,1,ntg,1,npg)
  call alloc_array2d(hwtpg_ex       ,1,ntg,1,npg)
end subroutine allocate_weight_midpoint_binary_excision
!
subroutine calc_weight_midpoint_binary_excision
  use grid_parameter,       only : nrg, ntg, npg 
  use weight_midpoint_grav, only : hwdrg, hwdtg, hwdpg, tzwdtg, tzwdpg
  use coordinate_grav_theta,   only : dthg
  use coordinate_grav_phi,     only : dphig, phig
  use grid_points_binary_excision, only : hphig_exin, hphig_exout
  use trigonometry_grav_theta, only : hsinthg
  implicit none
  integer :: irg, itg, ipg
  real(long) :: wgr, wgt, wgp, phi_in, phi_out, small = 1.0d-14
!
! --  weight of integlation for volume integral in GR coordinate, 
! --  with binary excision
!
  do itg = 1, ntg
    do irg = 1, nrg
      do ipg = 1, npg
!
        wgr = hwdrg(irg)
        wgt = hwdtg(itg)
        wgp = hwdpg(ipg)
        phi_in  = hphig_exin(irg,itg)
        phi_out = hphig_exout(irg,itg)
        if (phi_in.ge.small) then
          if (phi_in.ge.phig(ipg-1).and.phi_in.le.phig(ipg)) then 
            wgp = phig(ipg) - phi_in
          end if
          if (phi_out.ge.phig(ipg-1).and.phi_out.le.phig(ipg)) then 
            wgp = phi_out - phig(ipg-1)
          end if
          if (phig(ipg).le.phi_in.or.phig(ipg-1).ge.phi_out) then 
            wgr = 0.0d0
            wgt = 0.0d0
            wgp = 0.0d0
          end if
        end if
!
        hwrtpg_ex(irg,itg,ipg) = wgr*wgt*wgp
!
      end do
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      wgt = hwdtg(itg)
      wgp = hwdpg(ipg)
      hwtpg_ex(itg,ipg) = wgt*wgp
    end do
  end do
!
end subroutine calc_weight_midpoint_binary_excision
!
subroutine calc_weight_midpoint_binary_excision_hybrid
  use grid_parameter,       only : nrg, ntg, npg, ntgeq, nrf
  use weight_midpoint_grav, only : hwdrg, hwdtg, hwdpg, tzwdtg, tzwdpg, hwrtpg
  use coordinate_grav_theta,   only : dthg
  use coordinate_grav_phi,     only : dphig, phig
  use grid_points_binary_excision, only : hphig_exin, hphig_exout, &
  &                                       irg_exin, irg_exout
  use trigonometry_grav_theta, only : hsinthg
  implicit none
  integer :: irg, itg, ipg, nrg_exin, nrg_exout
  real(long) :: wgr, wgt, wgp, phi_in, phi_out, small = 1.0d-14
!
! --  weight of integlation for volume integral in GR coordinate, 
! --  with binary excision
!
  nrg_exin  = irg_exin(ntgeq,0) !- nrf/2
  nrg_exout = irg_exin(ntgeq,0) !+ nrf/2
  do itg = 1, ntg
    do irg = 1, nrg
      do ipg = 1, npg
!
        wgr = hwdrg(irg)
!        wgt = hwdtg(itg)
        wgt = hsinthg(itg)*dthg
!        if (irg.ge.nrg_exin.and.irg.le.nrg_exout) wgt = hsinthg(itg)*dthg
        wgp = hwdpg(ipg)
        phi_in  = hphig_exin(irg,itg)
        phi_out = hphig_exout(irg,itg)
        if (phi_in.ge.small) then
          if (phi_in.ge.phig(ipg-1).and.phi_in.le.phig(ipg)) then 
            wgp = phig(ipg) - phi_in
          end if
          if (phi_out.ge.phig(ipg-1).and.phi_out.le.phig(ipg)) then 
            wgp = phi_out - phig(ipg-1)
          end if
          if (phig(ipg).le.phi_in.or.phig(ipg-1).ge.phi_out) then 
            wgr = 0.0d0
            wgt = 0.0d0
            wgp = 0.0d0
          end if
        end if
!
        hwrtpg_ex(irg,itg,ipg) = wgr*wgt*wgp
!
      end do
    end do
  end do
!
  do ipg = 1, npg
    do itg = 1, ntg
      wgt = hwdtg(itg)
      wgp = hwdpg(ipg)
!!      hwtpg_ex(itg,ipg) = wgt*wgp
!!      hwtpgsf(it,ip) = hsinthg(it)*dthg*hwdpg(ip)
      do irg = 1, nrg
        wgt = hsinthg(itg)*dthg
        wgr = hwdrg(irg)
        hwrtpg(irg,itg,ipg) = wgr*wgt*wgp
      end do
    end do
  end do
!
end subroutine calc_weight_midpoint_binary_excision_hybrid
end module weight_midpoint_binary_excision
