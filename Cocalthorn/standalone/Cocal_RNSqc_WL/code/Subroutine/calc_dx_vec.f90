subroutine calc_dx_vec(cline,ir,it,ip,dxv)
  use phys_constant, only :  long
  use coordinate_grav_r, only  :   rg
  use def_matter, only : rs
  use trigonometry_grav_theta, only : hsinthg, hcosthg, sinthg, costhg
  use trigonometry_grav_phi, only : hsinphig, hcosphig, sinphig, cosphig
  use coordinate_grav_theta, only : dthginv
  use coordinate_grav_phi, only : dphiginv
  implicit none
  character(len=2), intent(in) :: cline
  integer,          intent(in) :: ir,it,ip
  real(long), pointer :: dxv(:)
  real(long) :: dprs, hprs, dtrs, htrs
!
  if (cline=="ph") then                    ! theta = constant, gridpoint
    dprs = rg(ir)*(rs(it,ip) - rs(it,ip-1))*dphiginv   ! at phi midpoint
    hprs = rg(ir)*(rs(it,ip) + rs(it,ip-1))*0.5d0      ! at phi midpoint

    dxv(1) = dprs*sinthg(it)*hcosphig(ip) - hprs*sinthg(it)*hsinphig(ip)
    dxv(2) = dprs*sinthg(it)*hsinphig(ip) + hprs*sinthg(it)*hcosphig(ip)
    dxv(3) = dprs*costhg(it)
  end if

  if (cline=="th") then                       ! phi = constant, gridpoint
    dtrs = rg(ir)*(rs(it,ip) - rs(it-1,ip))*dthginv   ! at theta midpoint
    htrs = rg(ir)*(rs(it,ip) + rs(it-1,ip))*0.5d0     ! at theta midpoint

    dxv(1) = dtrs*hsinthg(it)*cosphig(ip) + htrs*hcosthg(it)*cosphig(ip)
    dxv(2) = dtrs*hsinthg(it)*sinphig(ip) + htrs*hcosthg(it)*sinphig(ip)
    dxv(3) = dtrs*hcosthg(it)             - htrs*hsinthg(it)
  end if
!
end subroutine calc_dx_vec
