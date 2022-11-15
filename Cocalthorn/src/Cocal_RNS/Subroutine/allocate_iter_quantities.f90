subroutine allocate_iter_quantities(niq)
  use phys_constant, only : long
  use def_iter_quantities
  use def_bh_parameter, only : ome_bh, spin_bh
  use def_quantities
  use make_array_1d
  implicit none
  integer, intent(in) :: niq
  integer :: i
!
  call alloc_array1d(msec_x_oold,1,niq)
  call alloc_array1d(msec_x_old ,1,niq)
  call alloc_array1d(msec_f_oold,1,niq)
  call alloc_array1d(msec_f_old ,1,niq)
!
  i=0
  i=i+1; msec_x_oold(i) = ome_bh
  i=i+1; msec_x_oold(i) = spin_bh
!
  i=0
  i=i+1; msec_f_oold(i) = 0.0d0
  i=i+1; msec_f_oold(i) = 0.0d0
!
end subroutine allocate_iter_quantities
