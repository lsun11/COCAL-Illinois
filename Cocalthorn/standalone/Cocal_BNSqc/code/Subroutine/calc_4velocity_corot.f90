subroutine calc_4velocity_corot
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_matter_parameter, only : ber
  use def_matter, only : emd, utf, uxf, uyf, uzf
  use def_matter_velocity, only :  vxu, vyu, vzu
  implicit none
  real(long) :: emdfc, rhofc, prefc, hhfc, enefc
  integer :: ir, it, ip
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
!
        emdfc = emd(ir,it,ip)
        if (emdfc <= 1.0d-15) emdfc = 1.0d-15
        call peos_q2hprho(emdfc, hhfc, prefc, rhofc, enefc)
        utf(ir,it,ip) = hhfc/ber
        uxf(ir,it,ip) = utf(ir,it,ip)*vxu(ir,it,ip)
        uyf(ir,it,ip) = utf(ir,it,ip)*vyu(ir,it,ip)
        uzf(ir,it,ip) = utf(ir,it,ip)*vzu(ir,it,ip)
!
      end do
    end do
  end do
end subroutine calc_4velocity_corot
