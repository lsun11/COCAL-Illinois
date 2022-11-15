subroutine interpo_gr2gr_4th_mpt(fnc,cfn,rc,thc,phic)
  use phys_constant, only : long, pi
  use grid_parameter_interpo, only : nrg_itp, ntg_itp, npg_itp, rgin_itp
  use grid_parameter_binary_excision_interpo, only : ex_radius_itp
  use coordinate_grav_extended_interpo
  implicit none
  real(long), pointer     :: fnc(:,:,:)
  real(long), intent(out) :: cfn
  real(long) ::  rc, thc, phic
  real(long) ::  r4(4), th4(4), phi4(4), fr4(4), ft4(4), fp4(4)
  integer :: irg, itg, ipg, irgex, itgex, ipgex
  integer :: ir0, it0, ip0, irg0 , itg0 , ipg0, ii, jj, kk
  real(long), external :: lagint_4th
!
! --- Interpolation to a grid point (rc,thc,phic)
! --- using 4th order Lagrange formula.
!
  cfn = 0.0d0
  if (rc.lt.rgin_itp) return
  if (rc.gt.ex_radius_itp) stop 'interpo_gr2gr_4th_mpt'
!
  do irg = 0, nrg_itp+1
    if (rc.lt.rgex_itp(irg).and.rc.ge.rgex_itp(irg-1)) &
    &  ir0 = min0(irg-2,nrg_itp-3)
  end do
  do itg = 0, ntg_itp+1
    if (thc.lt.thgex_itp(itg).and.thc.ge.thgex_itp(itg-1)) it0 = itg-2
  end do
  do ipg = 0, npg_itp+1
    if (phic.lt.phigex_itp(ipg).and.phic.ge.phigex_itp(ipg-1)) ip0 = ipg-2
  end do
!
  do ii = 1, 4
    irg0 = ir0 + ii - 1
    itg0 = it0 + ii - 1
    ipg0 = ip0 + ii - 1
    r4(ii) = rgex_itp(irg0)
    th4(ii) = thgex_itp(itg0)
    phi4(ii) = phigex_itp(ipg0)
  end do
!
  do kk = 1, 4
    ipg0 = ip0 + kk - 1
    do jj = 1, 4
      itg0 = it0 + jj - 1
      do ii = 1, 4
        irg0 = ir0 + ii - 1
        irgex = irgex_r_itp(irg0)
        itgex = itgex_r_itp(itgex_th_itp(itg0),irg0)
        ipgex = ipgex_r_itp(ipgex_th_itp(ipgex_phi_itp(ipg0),itg0),irg0)
        fr4(ii) = fnc(irgex,itgex,ipgex)
      end do
      ft4(jj) = lagint_4th(r4,fr4,rc)
    end do
    fp4(kk) = lagint_4th(th4,ft4,thc)
  end do
  cfn = lagint_4th(phi4,fp4,phic)
!
end subroutine interpo_gr2gr_4th_mpt
