subroutine calc_matter_rhof
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_matter, only : emd, rhof
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
        rhof(ir,it,ip) = rhofc
!
      end do
    end do
  end do
end subroutine calc_matter_rhof
