subroutine plane_surf_int_fluid(cline, ia, ib, souf, surf_int)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   ntf, npf, nrf
  use coordinate_grav_r, only : drg, hrg
  use weight_midpoint_fluid, only : hwdpf, hwdtf, dthg
  use def_matter, only : rs
  implicit none
  character(len=2),  intent(in)  :: cline
  integer,           intent(in)  :: ia, ib
  real(long), pointer     :: souf(:,:)
  real(long), intent(out) :: surf_int
  real(long) :: hsou, hsurf
  integer    :: ipf, itf, irf, ir, it
!
  if (cline=="ph") then
    itf = ib
    surf_int = 0.0d0
    do irf = 1, ia
      do ipf = 1, npf
        hsou = souf(irf,ipf)
        hsurf = 0.5d0*(rs(itf,ipf) + rs(itf,ipf-1))
        surf_int = surf_int + hsou * hsurf * (hrg(irf)*drg(irf)*hwdpf(ipf))
      end do
    end do
  end if
!
  if (cline=="th") then
    surf_int = 0.0d0
    do irf = 1, ia
      ipf = ib
      do itf = 1, ntf
        hsou = souf(irf,itf)
        hsurf = 0.5d0*(rs(itf,ipf) + rs(itf-1,ipf))
        surf_int = surf_int + hsou * hsurf * (hrg(irf)*drg(irf)*dthg)
      end do

      ipf = ib + npf/2             ! ipf gridpoint
      do it = ntf+1, 2*ntf
        itf = ntf - (it - ntf) + 1   ! itf midpoint  ntf -> 1

        hsou = souf(irf,it)     ! Notice it not itf

        hsurf = 0.5d0*(rs(itf,ipf) + rs(itf-1,ipf))
        surf_int = surf_int + hsou * hsurf * (hrg(irf)*drg(irf)*dthg)
      end do
    end do
  end if
!
end subroutine plane_surf_int_fluid
