subroutine copy_def_peos_parameter_to_mpt(impt)
  use def_peos_parameter      !abc,abi,rhoi,qi,hi,nphase,rhoini_cgs,emdini_gcm1
  use def_peos_parameter_mpt
  use copy_array_static_1dto2d_mpt
  implicit none
  integer :: i, impt
!  
  call copy_arraystatic_1dto2d_mpt(impt, enei,   enei_,   0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, prei,   prei_,   0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, rhoi,   rhoi_,   0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, qi,     qi_,     0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, hi,     hi_,     0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, rhocgs, rhocgs_, 0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, enecgs, enecgs_, 0, nnteos)
  call copy_arraystatic_1dto2d_mpt(impt, precgs, precgs_, 0, nnteos)
!
  i=0
  i=i+1; def_peos_param_real_(i,impt) = rhoini_cgs
  i=i+1; def_peos_param_real_(i,impt) = rhoini_gcm1
  i=i+1; def_peos_param_real_(i,impt) = emdini_gcm1
  i=i+1; def_peos_param_real_(i,impt) = kappa_crust
  i=i+1; def_peos_param_real_(i,impt) = gamma_crust
!
  i=0
  i=i+1; def_peos_param_int_(i,impt) = nphase
!
end subroutine copy_def_peos_parameter_to_mpt
