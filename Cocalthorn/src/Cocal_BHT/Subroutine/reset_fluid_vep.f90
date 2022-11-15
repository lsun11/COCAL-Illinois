subroutine reset_fluid_vep
  use def_matter, only : vep, vepxf, vepyf, vepzf
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfeq, npfxzp, npfxzm, NS_shape
  use def_matter_parameter, only : ome
  use def_vector_phi, only : vec_phif
  implicit none
  integer :: ir, it, ip
!
! Impose symmetry of the shape to velocity potential
! NS_shape = JB : tri-axial configuration
! NS_shape = ML : axi-symmetric configuration
!
!  if (NS_shape.eq.'IB') then
! Impose xz and equatorial symmetry
    vep(0:nrf,0,0:npf) = 0.0d0
    vep(0:nrf,0:ntf,npfxzp) = 0.0d0
    vep(0:nrf,0:ntf,npfxzm) = 0.0d0
    vep(0:nrf,0:ntf,npf) = 0.0d0

    vep(0,0:ntf,0:npf) = 0.0d0

!   1/4 space:  y,z>0                          
    do it = 0, ntfeq 
      do ip = 0, npfxzm
        vep(0:nrf,ntf-it,    ip) =   vep(0:nrf,it,ip)   ! 1/4 space: y>0  z<0
        vep(0:nrf,    it,npf-ip) = - vep(0:nrf,it,ip)   ! 1/4 space: y<0  z>0
        vep(0:nrf,ntf-it,npf-ip) = - vep(0:nrf,it,ip)   ! 1/4 space: y<0  z<0
      end do
    end do
!  end if
!
end subroutine reset_fluid_vep
