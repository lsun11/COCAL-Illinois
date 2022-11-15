!  weight for numerical integration using mid-point rule
!______________________________________________
module weight_midpoint_fluid
  use phys_constant,           only : nnrf, nntf, nnpf, long
  use grid_parameter,          only : nrf, ntf, npf, nrg, ntg, npg
  use coordinate_grav_r,       only : drg, rg, hrg
  use coordinate_grav_theta,   only : dthg
  use coordinate_grav_phi,     only : dphig
  use trigonometry_grav_theta, only : hsinthg, sinthg
  use make_array_3d
  use make_array_2d
  implicit none
! weight
  Real(long)          ::  hwdrf(nnrf), hwdtf(nntf), hwdpf(nnpf)
  Real(long)          ::  tzwdrf(0:nnrf), siwdrf(0:nnrf)
  Real(long)          ::  siwdtf(0:nntf), tzwdpf(0:nnpf)
  Real(long)          ::  wdxf(0:nnrf)
  Real(long), pointer ::  hwrtpf(:,:,:), tzwrtpf(:,:,:), hwtpfsf(:,:)
  Real(long), pointer ::  siwrtpf(:,:,:), rtsiwrtpf(:,:,:)
contains
subroutine allocate_weight_midpoint_fluid
  implicit none
  call alloc_array2d(hwtpfsf,1,ntf,1,npf)
  call alloc_array3d(hwrtpf,1,nrf,1,ntf,1,npf)
  call alloc_array3d(tzwrtpf,0,nrf,1,ntf,1,npf)
  call alloc_array3d(siwrtpf,0,nrf,1,ntf,1,npf)
  call alloc_array3d(rtsiwrtpf,0,nrf,0,ntf,0,npf)
end subroutine allocate_weight_midpoint_fluid
subroutine weight_calc_midpoint_fluid
  implicit none
  integer             ::  ir, it, ip
  hwdrf(1:nrf) = hrg(1:nrf)**2*drg(1:nrf)
  hwdtf(1:ntf) = hsinthg(1:ntf)*dthg
  hwdpf(1:npf) = dphig
! Trapezoidal in r
  tzwdrf(1:nrf-1) = rg(1:nrf-1)**2*0.5d0*(drg(1:nrf-1) + drg(2:nrf))
  tzwdrf(0)   = 0.5d0*rg(0)**2*drg(1)
  tzwdrf(nrf) = 0.5d0*rg(nrf)**2*drg(nrf)
! Trapezoidal in r no Jacobian
  wdxf(1:nrf-1) = 0.5d0*(drg(1:nrf-1) + drg(2:nrf))
  wdxf(0)   = 0.5d0*drg(1)
  wdxf(nrf) = 0.5d0*drg(nrf)
! Simpson in r
  siwdrf(0) = 0.0d0
  do ir = 0, nrf-2, 2
    siwdrf(ir)   = siwdrf(ir) &
    &            + rg(ir  )**2*(2.0d0*drg(ir+1)-drg(ir+2)) &
    &                     *(drg(ir+1)+drg(ir+2))/(6.0d0*drg(ir+1))
    siwdrf(ir+1) = rg(ir+1)**2*(drg(ir+1)+drg(ir+2))**3 &
    &                   /(6.0d0*drg(ir+1)*drg(ir+2))
    siwdrf(ir+2) = rg(ir+2)**2*(2.0d0*drg(ir+2)-drg(ir+1)) &
    &                     *(drg(ir+1)+drg(ir+2))/(6.0d0*drg(ir+2))
  end do
! Simpson in theta
  siwdtf(0)   = 1.0d0/3.0d0*sinthg(0)*dthg
  siwdtf(1:ntf-1:2) = 4.0d0/3.0d0*sinthg(1:ntf-1:2)*dthg
  siwdtf(2:ntf-2:2) = 2.0d0/3.0d0*sinthg(2:ntf-2:2)*dthg
  siwdtf(ntf) = 1.0d0/3.0d0*sinthg(ntf)*dthg
! Trapezoidal in phi
  tzwdpf(1:npf-1) = dphig
  tzwdpf(0)   = 0.5d0*dphig 
  tzwdpf(npf) = 0.5d0*dphig
!
  do ip = 1, npf
    do it = 1, ntf
      hwtpfsf(it,ip) = hwdtf(it)*hwdpf(ip)
      tzwrtpf(0,it,ip) = tzwdrf(0)*hwdtf(it)*hwdpf(ip)
      siwrtpf(0,it,ip) = siwdrf(0)*hwdtf(it)*hwdpf(ip)
      do ir = 1, nrf
        hwrtpf(ir,it,ip) = hwdrf(ir)*hwdtf(it)*hwdpf(ip)
        tzwrtpf(ir,it,ip) = tzwdrf(ir)*hwdtf(it)*hwdpf(ip)
        siwrtpf(ir,it,ip) = siwdrf(ir)*hwdtf(it)*hwdpf(ip)
      end do
    end do
  end do
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
        rtsiwrtpf(ir,it,ip) = siwdrf(ir)*siwdtf(it)*tzwdpf(ip)
      end do
    end do
  end do
!
end subroutine weight_calc_midpoint_fluid
end module weight_midpoint_fluid
