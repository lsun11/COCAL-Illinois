!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_grav
  use phys_constant,           only : nnrg, nntg, nnpg, long
  use grid_parameter,          only : nrg, ntg, npg 
  use coordinate_grav_r,       only : drg, rg, hrg
  use coordinate_grav_theta,   only : dthg
  use coordinate_grav_phi,     only : dphig
  use trigonometry_grav_theta, only : sinthg, hsinthg
  use make_array_2d
  use make_array_3d
  implicit none
! weight
  Real(long)          ::  hwdrg(nnrg), hwdtg(nntg), hwdpg(nnpg)
  Real(long)          ::  wdtg(nntg)
  Real(long)          ::  tzwdrg(0:nnrg), tzwdtg(0:nntg), tzwdpg(0:nnpg)
  Real(long)          ::  wdxg(0:nnrg)
  Real(long),pointer  ::  hwrtpg(:,:,:), tzwrtpg(:,:,:), hwtpgsf(:,:)
contains
subroutine allocate_weight_midpoint_grav
  implicit none
  call alloc_array2d(hwtpgsf,1,ntg,1,npg)
  call alloc_array3d(hwrtpg,1,nrg,1,ntg,1,npg)
  call alloc_array3d(tzwrtpg,0,nrg,1,ntg,1,npg)
end subroutine allocate_weight_midpoint_grav
subroutine weight_calc_midpoint_grav
  implicit none
  integer             ::  ir, it, ip
  hwdrg(1:nrg) = hrg(1:nrg)**2*drg(1:nrg)
  hwdtg(1:ntg) = hsinthg(1:ntg)*dthg
  wdtg(1:ntg)  = dthg
  hwdpg(1:npg) = dphig
! Trapezoidal in r, theta, phi
  tzwdrg(1:nrg-1) = rg(1:nrg-1)**2*0.5d0*(drg(1:nrg-1) + drg(2:nrg))
  tzwdrg(0)   = 0.5d0*rg(0)**2*drg(1)
  tzwdrg(nrg) = 0.5d0*rg(nrg)**2*drg(nrg)
  tzwdtg(1:ntg-1) = sinthg(1:ntg-1)*dthg
  tzwdtg(0)   = 0.5d0*sinthg(0)*dthg
  tzwdtg(ntg) = 0.5d0*sinthg(ntg)*dthg
  tzwdpg(1:npg-1) = dphig
  tzwdpg(0)   = 0.5d0*dphig
  tzwdpg(npg) = 0.5d0*dphig
! Trapezoidal in r no Jacobian
  wdxg(1:nrg-1) = 0.5d0*(drg(1:nrg-1) + drg(2:nrg))
  wdxg(0)   = 0.5d0*drg(1)
  wdxg(nrg) = 0.5d0*drg(nrg)
!
  do ip = 1, npg
    do it = 1, ntg
      hwtpgsf(it,ip) = hwdtg(it)*hwdpg(ip)
      tzwrtpg(0,it,ip) = tzwdrg(0)*hwdtg(it)*hwdpg(ip)
      do ir = 1, nrg
        hwrtpg(ir,it,ip) = hwdrg(ir)*hwdtg(it)*hwdpg(ip)
        tzwrtpg(ir,it,ip) = tzwdrg(ir)*hwdtg(it)*hwdpg(ip)
      end do
    end do
  end do
!
end subroutine weight_calc_midpoint_grav
subroutine weight_calc_midpoint_grav_th4th
  implicit none
  integer             ::  ir, it, ip
!
  hwdrg(1:nrg) = hrg(1:nrg)**2*drg(1:nrg)
  do it = 1, ntg, 4
    hwdtg(it  ) = hsinthg(it  )*dthg*13.0d0/12.0d0
    hwdtg(it+1) = hsinthg(it+1)*dthg*11.0d0/12.0d0
    hwdtg(it+2) = hsinthg(it+2)*dthg*11.0d0/12.0d0
    hwdtg(it+3) = hsinthg(it+3)*dthg*13.0d0/12.0d0

    wdtg(it  ) = dthg*13.0d0/12.0d0
    wdtg(it+1) = dthg*11.0d0/12.0d0
    wdtg(it+2) = dthg*11.0d0/12.0d0
    wdtg(it+3) = dthg*13.0d0/12.0d0
  end do
!!  do it = 1, ntg, 3
!!    hwdtg(it  ) = hsinthg(it  )*dthg*9.0d0/8.0d0
!!    hwdtg(it+1) = hsinthg(it+1)*dthg*3.0d0/4.0d0
!!    hwdtg(it+2) = hsinthg(it+2)*dthg*9.0d0/8.0d0
!!  end do
  hwdpg(1:npg) = dphig
! Trapezoidal in r, theta, phi
  tzwdrg(1:nrg-1) = rg(1:nrg-1)**2*0.5d0*(drg(1:nrg-1) + drg(2:nrg))
  tzwdrg(0)   = 0.5d0*rg(0)**2*drg(1)
  tzwdrg(nrg) = 0.5d0*rg(nrg)**2*drg(nrg)
  tzwdtg(1:ntg-1) = sinthg(1:ntg-1)*dthg
  tzwdtg(0)   = 0.5d0*sinthg(0)*dthg
  tzwdtg(ntg) = 0.5d0*sinthg(ntg)*dthg
  tzwdpg(1:npg-1) = dphig
  tzwdpg(0)   = 0.5d0*dphig
  tzwdpg(npg) = 0.5d0*dphig
! Trapezoidal in r no Jacobian
  wdxg(1:nrg-1) = 0.5d0*(drg(1:nrg-1) + drg(2:nrg))
  wdxg(0)   = 0.5d0*drg(1)
  wdxg(nrg) = 0.5d0*drg(nrg)
!
  do ip = 1, npg
    do it = 1, ntg
      hwtpgsf(it,ip) = hwdtg(it)*hwdpg(ip)
      tzwrtpg(0,it,ip) = tzwdrg(0)*hwdtg(it)*hwdpg(ip)
      do ir = 1, nrg
        hwrtpg(ir,it,ip) = hwdrg(ir)*hwdtg(it)*hwdpg(ip)
        tzwrtpg(ir,it,ip) = tzwdrg(ir)*hwdtg(it)*hwdpg(ip)
      end do
    end do
  end do
!
end subroutine weight_calc_midpoint_grav_th4th
end module weight_midpoint_grav
