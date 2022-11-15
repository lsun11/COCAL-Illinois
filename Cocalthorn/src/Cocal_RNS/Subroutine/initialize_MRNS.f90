subroutine initialize_MRNS
  use phys_constant, only  : long
  use grid_parameter, only : nrf, ntf, npf
  use def_matter, only : emd, rs, utf, uxf, uyf, uzf
  use def_matter_parameter, only : ome, ber
  use def_vector_phi, only : vec_phif
  implicit none
  real(long) :: hh, pre, rho, ene, qq
  integer :: irf, itf, ipf
!
! --- initialize matter variables
!
  call calc_vector_phi_matter(1)
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf
        qq = emd(irf,itf,ipf)
        call peos_q2hprho(qq, hh, pre, rho, ene)
        utf(irf,itf,ipf) = hh/ber
        uxf(irf,itf,ipf) = vec_phif(irf,itf,ipf,1)*ome
        uyf(irf,itf,ipf) = vec_phif(irf,itf,ipf,2)*ome
        uzf(irf,itf,ipf) = 0.0d0
      end do
    end do
  end do
!
end subroutine initialize_MRNS
