subroutine copy_def_bh_parameter_from_mpt(impt)
  use def_bh_parameter_mpt
  use def_bh_parameter
  implicit none
  integer :: i, impt
!  
  i=0
  i=i+1; ome_bh  = def_bh_param_real_(i,impt)
  i=i+1; spin_bh = def_bh_param_real_(i,impt)
  i=i+1; alph_bh = def_bh_param_real_(i,impt)
  i=i+1; psi_bh  = def_bh_param_real_(i,impt)
  i=i+1; th_spin_bh_deg  =  def_bh_param_real_(i,impt)
  i=i+1; phi_spin_bh_deg =  def_bh_param_real_(i,impt)
  i=i+1; th_spin_bh  =  def_bh_param_real_(i,impt)
  i=i+1; phi_spin_bh =  def_bh_param_real_(i,impt)
  i=i+1; rgin_deform =  def_bh_param_real_(i,impt)
!
  i=i+1; mass_pBH    = def_bh_param_real_(i,impt)
  i=i+1; mom_pBH(1)  = def_bh_param_real_(i,impt)
  i=i+1; mom_pBH(2)  = def_bh_param_real_(i,impt)
  i=i+1; mom_pBH(3)  = def_bh_param_real_(i,impt)
  i=i+1; spin_pBH(1) = def_bh_param_real_(i,impt)
  i=i+1; spin_pBH(2) = def_bh_param_real_(i,impt)
  i=i+1; spin_pBH(3) = def_bh_param_real_(i,impt)
!
  i=0
  i=i+1; bh_bctype  = def_bh_param_char2_(i,impt)
  i=i+1; bh_sptype  = def_bh_param_char2_(i,impt)
  i=i+1; bh_soltype = def_bh_param_char2_(i,impt)
!
end subroutine copy_def_bh_parameter_from_mpt
