module radial_green_fn_hret_mi_hadv
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  bsjy_hrmiha(:,:,:,:)
  real(long), pointer  ::  sbsjy_hrmiha(:,:,:), sbsjyp_hrmiha(:,:,:)
!
  contains
  subroutine allocate_radial_green_fn_hret_mi_hadv
    use grid_parameter, only : nrg, nlg
    use make_array_3d
    use make_array_4d
    implicit none 
    call alloc_array4d(bsjy_hrmiha,1,nrg,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjy_hrmiha,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjyp_hrmiha,0,nlg,0,nlg,0,nrg)
  end subroutine allocate_radial_green_fn_hret_mi_hadv
end module radial_green_fn_hret_mi_hadv
