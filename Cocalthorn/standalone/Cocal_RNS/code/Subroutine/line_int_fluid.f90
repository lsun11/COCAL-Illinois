subroutine line_int_fluid(cline, ia, ib, souf, line_int)
  use phys_constant, only  :   long, pi
  use grid_parameter, only  :   ntf, npf
  use weight_midpoint_fluid, only : hwdpf, hwdtf, dthg
  use interface_calc_dx_vec
  use make_array_1d
  implicit none
  character(len=2),  intent(in)  :: cline
  integer,           intent(in)  :: ia, ib
  real(long), pointer     :: souf(:,:)
  real(long), intent(out) :: line_int
  real(long) :: hsou
  integer    :: ipf, itf, irf, it, jj
  real(long), pointer     :: dxv(:)

  call alloc_array1d(dxv, 1, 3)

  if (cline=="ph") then
    irf = ia
    itf = ib
    line_int = 0.0d0
    do ipf = 1, npf
      call calc_dx_vec(cline, irf, itf, ipf, dxv)
      hsou = souf(ipf,1)*dxv(1) + souf(ipf,2)*dxv(2) + souf(ipf,3)*dxv(3)
      line_int = line_int + hsou * hwdpf(ipf)
    end do
  end if
!
  if (cline=="th") then
    irf = ia
    ipf = ib
    line_int = 0.0d0
    do itf = 1, ntf
      call calc_dx_vec(cline, irf, itf, ipf, dxv)
!      if (itf==(ntf-2))  write(6,'(1p,3e23.15)')   (dxv(jj), jj=1,3)

      hsou = souf(itf,1)*dxv(1) + souf(itf,2)*dxv(2) + souf(itf,3)*dxv(3)
      line_int = line_int + hsou * dthg
    end do
!    write(6,'(a34,1p,e23.15)')  "Only for theta [0,pi] => line_int=", line_int

    do it = ntf+1, 2*ntf
      itf = ntf - (it - ntf) + 1   ! itf midpoint
      ipf = ib + npf/2             ! ipf gridpoint
      call calc_dx_vec(cline, irf, itf, ipf, dxv)
      dxv(1:3) = -dxv(1:3)
!      if (itf==(ntf-2))  write(6,'(1p,3e23.15)')   (dxv(jj), jj=1,3)

      hsou = souf(it,1)*dxv(1) + souf(it,2)*dxv(2) + souf(it,3)*dxv(3)
      line_int = line_int + hsou * dthg
    end do
  end if
!
  deallocate(dxv)

end subroutine line_int_fluid
