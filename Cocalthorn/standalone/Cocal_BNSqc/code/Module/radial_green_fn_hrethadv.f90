module radial_green_fn_hrethadv
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  bsjy_hrha(:,:,:,:)
  real(long), pointer  ::  sbsjy_hrha(:,:,:), sbsjyp_hrha(:,:,:)
!
  contains
  subroutine allocate_radial_green_fn_hrethadv
    use grid_parameter, only : nrg, nlg
    use make_array_3d
    use make_array_4d
    implicit none 
    call alloc_array4d(bsjy_hrha,1,nrg,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjy_hrha,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjyp_hrha,0,nlg,0,nlg,0,nrg)
  end subroutine allocate_radial_green_fn_hrethadv
end module radial_green_fn_hrethadv
