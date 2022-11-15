subroutine copy_def_matter_parameter_from_mpt(impt)
  use def_matter_parameter_mpt
  use def_matter_parameter
  implicit none
  integer :: i, impt
!  
  i=0
  i=i+1; pinx      = def_matter_param_real_(i,impt)  ! 1 
  i=i+1; emdc      = def_matter_param_real_(i,impt)  ! 2
  i=i+1; ome       = def_matter_param_real_(i,impt)  ! 3
  i=i+1; ber       = def_matter_param_real_(i,impt)  ! 4
  i=i+1; radi      = def_matter_param_real_(i,impt)  ! 5
  i=i+1; A2DR      = def_matter_param_real_(i,impt)  ! 6
  i=i+1; DRAT_A2DR = def_matter_param_real_(i,impt)  ! 7
  i=i+1; index_DR  = def_matter_param_real_(i,impt)  ! 8
  i=i+1; index_DRq = def_matter_param_real_(i,impt)  ! 9
  i=i+1; B2DR      = def_matter_param_real_(i,impt)  ! 10
  i=i+1; DRAT_B2DR = def_matter_param_real_(i,impt)  ! 11
  i=i+1; index_DRp = def_matter_param_real_(i,impt)  ! 12
  i=i+1; omespx    = def_matter_param_real_(i,impt)  ! 13
  i=i+1; omespy    = def_matter_param_real_(i,impt)  ! 14
  i=i+1; omespz    = def_matter_param_real_(i,impt)  ! 15
  i=i+1; confpow   = def_matter_param_real_(i,impt)  ! 16
  i=i+1; velx      = def_matter_param_real_(i,impt)  ! 17
  i=i+1; delome    = def_matter_param_real_(i,impt)  ! 18
  i=i+1; delvel    = def_matter_param_real_(i,impt)  ! 19
!
  i=0
  i=i+1; ROT_LAW = def_matter_param_char2_(i,impt) 
!
end subroutine copy_def_matter_parameter_from_mpt
