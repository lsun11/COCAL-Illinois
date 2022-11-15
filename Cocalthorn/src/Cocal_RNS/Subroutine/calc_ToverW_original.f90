subroutine calc_ToverW
  use phys_constant, only  :   long
  use def_matter_parameter, only  :   ome, radi
  use def_quantities
  implicit none
!
  omega = ome/radi
  T_kinene = 0.5d0*omega*angmom
  W_gravene = propermass + T_kinene - admmass
  ToverW = T_kinene/dabs(W_gravene)
  I_inertia = angmom/omega
!
end subroutine calc_ToverW
