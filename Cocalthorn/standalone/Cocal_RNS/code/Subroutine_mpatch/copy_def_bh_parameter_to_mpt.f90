subroutine copy_def_bh_parameter_to_mpt(impt)
  use def_bh_parameter_mpt
  use def_bh_parameter
  implicit none
  integer :: i, impt
!  
  i=0
  i=i+1; def_bh_param_real_(i,impt) = ome_bh
  i=i+1; def_bh_param_real_(i,impt) = spin_bh
  i=i+1; def_bh_param_real_(i,impt) = alph_bh
  i=i+1; def_bh_param_real_(i,impt) = psi_bh
  i=i+1; def_bh_param_real_(i,impt) = th_spin_bh_deg
  i=i+1; def_bh_param_real_(i,impt) = phi_spin_bh_deg
  i=i+1; def_bh_param_real_(i,impt) = th_spin_bh
  i=i+1; def_bh_param_real_(i,impt) = phi_spin_bh
  i=i+1; def_bh_param_real_(i,impt) = rgin_deform
!
  i=i+1; def_bh_param_real_(i,impt) = mass_pBH
  i=i+1; def_bh_param_real_(i,impt) = mom_pBH(1)
  i=i+1; def_bh_param_real_(i,impt) = mom_pBH(2)
  i=i+1; def_bh_param_real_(i,impt) = mom_pBH(3)
  i=i+1; def_bh_param_real_(i,impt) = spin_pBH(1)
  i=i+1; def_bh_param_real_(i,impt) = spin_pBH(2)
  i=i+1; def_bh_param_real_(i,impt) = spin_pBH(3)
!
  i=0
  i=i+1; def_bh_param_char2_(i,impt) =  bh_bctype
  i=i+1; def_bh_param_char2_(i,impt) =  bh_sptype
  i=i+1; def_bh_param_char2_(i,impt) =  bh_soltype
!
end subroutine copy_def_bh_parameter_to_mpt
