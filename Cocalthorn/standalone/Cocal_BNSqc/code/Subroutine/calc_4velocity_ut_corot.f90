subroutine calc_4velocity_ut_corot
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_matter_parameter, only : ber
  use def_matter, only : emd, utf, utg
  use interface_interpo_fl2gr
  implicit none
  real(long) :: emdfc, rhofc, prefc, hhfc, ene
  integer :: ir, it, ip
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
!
        emdfc = emd(ir,it,ip)
        if (emdfc <= 1.0d-15) emdfc = 1.0d-15
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, ene)
        utf(ir,it,ip) = hhfc/ber
!
      end do
    end do
  end do
  call interpo_fl2gr(utf,utg)
end subroutine calc_4velocity_ut_corot
