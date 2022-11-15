module def_dvphi_mpt
  use phys_constant, only : nnmpt
  implicit none
  real(8) :: dphiu_(3,3,nnmpt)
! --  dphiu(a,b) = D_b\phi^a
! --  dphiu(1,2) = -1.0, dphiu(2,1) = 1.0 
end module def_dvphi_mpt
