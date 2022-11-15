module radial_green_fn_hrethadv_homosol
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  sbsjy_hrha_ho(:,:,:), sbsjyp_hrha_ho(:,:,:)
!
  contains
  subroutine allocate_radial_green_fn_hrethadv_homosol
    use grid_parameter, only : nrg, nlg
    use make_array_3d
    implicit none 
    call alloc_array3d(sbsjy_hrha_ho,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjyp_hrha_ho,0,nlg,0,nlg,0,nrg)
  end subroutine allocate_radial_green_fn_hrethadv_homosol
end module radial_green_fn_hrethadv_homosol
