subroutine printout_iterqt_BNS_mpt(x_vector, f_vector, im,id,is)
  use phys_constant,  only : nmpt, g, c, solmas, long
  use def_quantities_mpt
  use def_matter_parameter_mpt
  use def_binary_parameter_mpt
  implicit none
  real(long), pointer :: x_vector(:), f_vector(:)
  real(8) :: MM = solmas, LL = g*solmas/c**2, TT = g*solmas/c**3
  integer :: im,id,is
!
  if(id.eq.1.and.im.eq.1)  then
    call length_conversion(2, surf_param_real_(2,1), bd)    !Target distance in km: surf_param_real_(2,1) 
    bsepa   = sepa_(1)
    write(6,'(a41,1p,e23.15)') "Binary separation in COCAL units        :", bsepa
    write(6,'(a41,1p,e23.15)') "Target binary separation in km          :", surf_param_real_(2,1)
    write(6,'(a41,1p,e23.15)') "Target binary separation in COCAL units :", bd
  end if
  
!
end subroutine printout_iterqt_BNS_mpt
