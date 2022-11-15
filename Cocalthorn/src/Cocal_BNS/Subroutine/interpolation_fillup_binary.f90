subroutine interpolation_fillup_binary(fnc)
  use phys_constant, only : long, pi
  use grid_parameter, only : ntg, npg, rgin
  use grid_parameter_binary_excision, only : ex_rgmid, ex_radius 
  use grid_points_binary_excision, only : rb, thb, phib, irg_exin, irg_exout
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use interface_interpo_gr2gr_4th
  implicit none
  real(long), pointer :: fnc(:,:,:)
  real(long) :: rc, thc, phic, cfn
  integer :: irg, itg, ipg, nrg_exin, nrg_exout
!
  do ipg = 0, npg
    do itg = 0, ntg
      nrg_exin  = irg_exin(itg,ipg)  + 1
      nrg_exout = irg_exout(itg,ipg) - 1
      if (nrg_exin.eq.0) cycle
      do irg = nrg_exin, nrg_exout
        rc = rb(irg,itg,ipg)
        thc = thb(irg,itg,ipg)
        phic = dmod(phib(irg,itg,ipg)+pi,2.0d0*pi)
        if (rc.lt.rgin) cycle
        if (rc.gt.ex_radius) stop ' interpolation_fillup_binary '
        call interpo_gr2gr_4th(fnc,cfn,rc,thc,phic)
        fnc(irg,itg,ipg) = cfn
      end do
    end do
  end do
!
end subroutine interpolation_fillup_binary
