module def_bh_parameter_mpt
  use phys_constant, only : long
  implicit none
  real(long), pointer :: def_bh_param_real_(:,:)
  character(len=2), pointer :: def_bh_param_char2_(:,:)
end module def_bh_parameter_mpt
