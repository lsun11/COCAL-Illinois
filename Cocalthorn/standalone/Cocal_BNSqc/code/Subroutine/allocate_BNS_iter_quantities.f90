subroutine allocate_BNS_iter_quantities(niq)
  use phys_constant, only : long
  use def_iter_quantities
  use make_array_1d
  implicit none
  integer, intent(in) :: niq
!
  call alloc_array1d(msec_x_oold,1,niq)
  call alloc_array1d(msec_x_old ,1,niq)
  call alloc_array1d(msec_f_oold,1,niq)
  call alloc_array1d(msec_f_old ,1,niq)
!
end subroutine allocate_BNS_iter_quantities
