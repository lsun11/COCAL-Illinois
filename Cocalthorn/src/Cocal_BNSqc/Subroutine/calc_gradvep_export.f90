subroutine calc_gradvep_export(potf,potxf,potyf,potzf,rs)
  use phys_constant, only  :   long
  use grid_parameter
!  use def_matter_parameter, only : ome
!  use def_vector_phi, only : vec_phif
!  use def_matter, only : vep
  use interface_flgrad_4th_gridpoint
  use interface_flgrad_2nd_gridpoint_export
  implicit none
  real(long), pointer :: potf(:,:,:), potxf(:,:,:), potyf(:,:,:), potzf(:,:,:), rs(:,:)
  real(long) :: dxvep, dyvep, dzvep
  integer    :: irf, itf, ipf
!
  do ipf = 0, npf
    do itf = 0, ntf
      do irf = 0, nrf-1
        call flgrad_2nd_gridpoint_export(potf,dxvep,dyvep,dzvep,irf,itf,ipf,rs)
        potxf(irf,itf,ipf) = dxvep
        potyf(irf,itf,ipf) = dyvep
        potzf(irf,itf,ipf) = dzvep
      end do
    end do
  end do
!
  do ipf = 0, npf
    do itf = 0, ntf
      potxf(nrf,itf,ipf) = 2.0d0*potxf(nrf-1,itf,ipf) - potxf(nrf-2,itf,ipf)
      potyf(nrf,itf,ipf) = 2.0d0*potyf(nrf-1,itf,ipf) - potyf(nrf-2,itf,ipf)
      potzf(nrf,itf,ipf) = 2.0d0*potzf(nrf-1,itf,ipf) - potzf(nrf-2,itf,ipf)
    end do
  end do

!  write(6,*) " "
!  write(6,*) "gradvep x component:", potxf(nrf,ntf/2,0), potxf(0,0,0), potxf(nrf,ntf/2,npf/2)
!  write(6,*) "gradvep y component:", potyf(nrf,ntf/2,0), potyf(0,0,0), potyf(nrf,ntf/2,npf/2)

end subroutine calc_gradvep_export
