subroutine calc_3velocity_corot
  use phys_constant, only : long
  use grid_parameter, only : nrf, ntf, npf
  use def_matter_parameter, only : ome
  use def_matter_velocity, only : vxu, vyu, vzu
  use def_vector_phi, only : vec_phif
  implicit none
  integer :: ir, it, ip
!
  do ip = 0, npf
    do it = 0, ntf
      do ir = 0, nrf
!
        vxu(ir,it,ip) = ome*vec_phif(ir,it,ip,1)
        vyu(ir,it,ip) = ome*vec_phif(ir,it,ip,2)
        vzu(ir,it,ip) = 0.0d0
!
      end do
    end do
  end do
end subroutine calc_3velocity_corot
