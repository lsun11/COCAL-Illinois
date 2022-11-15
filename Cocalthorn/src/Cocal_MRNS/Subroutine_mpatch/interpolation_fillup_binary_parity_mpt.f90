subroutine interpolation_fillup_binary_parity_mpt(fnc,fnc_ex,parchar)
  use phys_constant, only : long, pi
  use grid_parameter, only : ntg, npg, rgin
  use grid_parameter_binary_excision, only : ex_rgmid, ex_radius 
  use grid_points_binary_excision, only : rb, thb, phib, irg_exin, irg_exout
  use coordinate_grav_r, only : rg
  use coordinate_grav_theta, only : thg
  use coordinate_grav_phi, only : phig
  use interface_interpo_gr2gr_4th_mpt
  implicit none
  real(long), pointer :: fnc(:,:,:), fnc_ex(:,:,:)
  real(long) :: rc, thc, phic, cfn
  real(long) :: pari
  integer :: irg, itg, ipg, nrg_exin, nrg_exout
  character(len=2) :: parchar
!
  if (parchar.eq.'ev') pari = + 1.0d0
  if (parchar.eq.'od') pari = - 1.0d0
!
! fnc is in a patch impt, fnc_ex is in a patch impt_ex
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
        call interpo_gr2gr_4th_mpt(fnc_ex,cfn,rc,thc,phic)
        fnc(irg,itg,ipg) = cfn*pari
      end do
    end do
  end do
!
end subroutine interpolation_fillup_binary_parity_mpt
