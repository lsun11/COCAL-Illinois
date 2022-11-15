module radial_green_fn_helmholtz
  use phys_constant, only : long
  implicit none
  real(long), pointer  ::  bsjy(:,:,:,:), sbsjy(:,:,:), sbsjyp(:,:,:)
!
  contains
  subroutine allocate_radial_green_fn_helmholtz
    use grid_parameter, only : nrg, nlg
    use make_array_3d
    use make_array_4d
    implicit none 
    call alloc_array4d(bsjy,1,nrg,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjy,0,nlg,0,nlg,0,nrg)
    call alloc_array3d(sbsjyp,0,nlg,0,nlg,0,nrg)
  end subroutine allocate_radial_green_fn_helmholtz
end module radial_green_fn_helmholtz
