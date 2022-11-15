subroutine adjust_matter_to_surface(new_rsf)
  use phys_constant, only : long, nmpt
  use def_matter, only : emd, rs
  use def_matter_parameter, only : emdc
  use grid_parameter, only : nrf, ntf, npf, &
  &                          ntfeq, ntfxy, npfyzp, npfxzp, npfxzm, &
  &                          ratio, NS_shape, EQ_point, r_surf
  implicit none
  integer :: it, ip
  real(long) :: new_rsf, rat


  rat = new_rsf/rs(ntfeq,npfxzp)

  rs(0:ntf,0:npf) = rat*rs(0:ntf,0:npf)                   

!
end subroutine adjust_matter_to_surface
