subroutine copy_def_matter_parameter_to_mpt(impt)
  use def_matter_parameter_mpt
  use def_matter_parameter
  implicit none
  integer :: i, impt
!  
  i=0
  i=i+1; def_matter_param_real_(i,impt) = pinx
  i=i+1; def_matter_param_real_(i,impt) = emdc
  i=i+1; def_matter_param_real_(i,impt) = ome
  i=i+1; def_matter_param_real_(i,impt) = ber
  i=i+1; def_matter_param_real_(i,impt) = radi
  i=i+1; def_matter_param_real_(i,impt) = A2DR
  i=i+1; def_matter_param_real_(i,impt) = DRAT_A2DR 
  i=i+1; def_matter_param_real_(i,impt) = index_DR
  i=i+1; def_matter_param_real_(i,impt) = index_DRq 
  i=i+1; def_matter_param_real_(i,impt) = B2DR
  i=i+1; def_matter_param_real_(i,impt) = DRAT_B2DR 
  i=i+1; def_matter_param_real_(i,impt) = index_DRp
  i=i+1; def_matter_param_real_(i,impt) = omespx
  i=i+1; def_matter_param_real_(i,impt) = omespy
  i=i+1; def_matter_param_real_(i,impt) = omespz
  i=i+1; def_matter_param_real_(i,impt) = confpow
  i=i+1; def_matter_param_real_(i,impt) = velx
  i=i+1; def_matter_param_real_(i,impt) = delome
  i=i+1; def_matter_param_real_(i,impt) = delvel 
!
  i=0
  i=i+1; def_matter_param_char2_(i,impt) = ROT_LAW
!
end subroutine copy_def_matter_parameter_to_mpt
