subroutine reset_kerr_parameter
  use phys_constant, only : long, pi
  use def_bh_parameter, only : ome_bh, mass_pBH, m_kerr, a_kerr
  implicit none
  m_kerr = mass_pBH
  a_kerr = ome_bh
end subroutine reset_kerr_parameter
