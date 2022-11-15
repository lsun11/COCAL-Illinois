subroutine allocate_matter_4velocity
  use phys_constant, only : long
  use grid_parameter
  use def_matter, only  : utf, uxf, uyf, uzf, utdf, uxdf, uydf, uzdf, &
  &                       utg, uxg, uyg, uzg, utdg, uxdg, uydg, uzdg
  use make_array_3d
  implicit none
!
  call alloc_array3d(utf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(uxf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(uyf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(uzf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(utdf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(uxdf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(uydf, 0,nrf, 0,ntf, 0,npf)
  call alloc_array3d(uzdf, 0,nrf, 0,ntf, 0,npf)
!
  call alloc_array3d(utg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(uxg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(uyg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(uzg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(utdg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(uxdg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(uydg, 0,nrg, 0,ntg, 0,npg)
  call alloc_array3d(uzdg, 0,nrg, 0,ntg, 0,npg)
!
end subroutine allocate_matter_4velocity
