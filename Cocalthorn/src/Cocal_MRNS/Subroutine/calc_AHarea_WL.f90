subroutine calc_AHarea_WL
  use phys_constant, only  : long, pi
  use grid_parameter, only : ntg, npg
  use def_quantities, only : app_hor_area_bh
  use interface_source_AHarea_WL
  use interface_surf_int_grav_rg_noweights
  use make_array_2d
  implicit none
  real(long) :: surf
  real(long), pointer :: sousg(:,:)
!
  call alloc_array2d(sousg,0,ntg,0,npg)
!
  call source_AHarea_WL(sousg)          ! contains the weights
  call surf_int_grav_rg_noweights(sousg,surf)
  app_hor_area_bh = surf

  write(6,'(a14,1p,e15.7)') '*** AH area = ', app_hor_area_bh

  deallocate(sousg)
end subroutine calc_AHarea_WL
