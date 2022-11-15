subroutine reset_fluid_gradvep
  use def_matter, only : vep, vepxf, vepyf, vepzf
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfeq, npfxzp, npfxzm, NS_shape
  use def_metric_on_SFC_CF, only : bvydf
  use def_matter_parameter, only : ome
  use def_vector_phi, only : vec_phif
  implicit none
  integer :: irf, itf, ipf
!
!
! origin
  vepxf(0,0:ntf,0:npf) = 0.0d0
!  vepyf(0,0:ntf,0:npf) = bvydf(0,0,0) + ome*vec_phif(0,0,0,2)
  vepzf(0,0:ntf,0:npf) = 0.0d0

! x-axis
  vepxf(0:nrf,ntfeq,0)      = 0.0d0
  vepxf(0:nrf,ntfeq,npf)    = 0.0d0
  vepxf(0:nrf,ntfeq,npfxzm) = 0.0d0
  vepzf(0:nrf,ntfeq,0)      = 0.0d0
  vepzf(0:nrf,ntfeq,npf)    = 0.0d0
  vepzf(0:nrf,ntfeq,npfxzm) = 0.0d0

! xy plane
  do irf = 1, nrf
    do ipf = 0, npf
      vepzf(irf,ntfeq,ipf) = 0.0d0
    end do
  end do


!
end subroutine reset_fluid_gradvep
