module integrability_fnc_MHD
  use phys_constant, only : long
  implicit none
  real(long) :: MHDfnc_PSI, MHDfnc_dPSI, MHDfnc_d2PSI
  real(long) :: MHDfnc_At, MHDfnc_dAt, MHDfnc_d2At
  real(long) :: MHDfnc_Lambda,    MHDfnc_dLambda
  real(long) :: MHDfnc_Lambda_GS, MHDfnc_dLambda_GS
  real(long) :: MHDfnc_Lambda_t
  real(long) :: MHDfnc_Lambda_phi, MHDfnc_dLambda_phi
  real(long) :: MHDpar_apsi, MHDidx_p, MHDpar_a, MHDidx_k
  real(long) :: MHDpar_charge, MHDidx_q
  real(long) :: MHDpar_Lc, MHDidx_s
  real(long) :: Aphi_max_vol, Aphi_max_surf
end module integrability_fnc_MHD
