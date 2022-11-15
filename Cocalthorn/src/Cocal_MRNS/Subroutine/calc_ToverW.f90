subroutine calc_ToverW
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_matter_parameter, only : ome, radi
  use def_quantities
  use def_quantities_derived, only : omega
  use make_array_3d
  use interface_source_mp_minus_madm
  use interface_vol_int_grav
  use interface_vol_int_fluid
  implicit none
  real(long) :: fac2pi
  real(long) :: volg, volf
  real(long), pointer :: soug(:,:,:), souf(:,:,:)
!
  call alloc_array3d(soug, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(souf, 0, nrf, 0, ntf, 0, npf)
!
  omega = ome/radi
  T_kinene = 0.5d0*omega*angmom
!!  W_gravene = propermass + T_kinene - admmass
  call source_mp_minus_madm(soug,souf)
  call vol_int_grav(soug,volg)
  call vol_int_fluid(souf,volf)
!
  fac2pi = 0.5d0/pi
  W_gravene = fac2pi*(radi*volg + radi**3*volf) + T_kinene
!
  ToverW = T_kinene/dabs(W_gravene)
  I_inertia = angmom/omega
!
  deallocate(soug)
  deallocate(souf)
end subroutine calc_ToverW
