subroutine calc_ToverW_peos
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_matter, only : omef
  use def_matter_parameter, only : ome, radi
  use def_quantities
  use def_quantities_derived, only : omega
  use make_array_3d
  use interface_source_mp_minus_madm_peos
  use interface_vol_int_grav
  use interface_vol_int_fluid
  use interface_source_ang_mom_peos
  implicit none
  real(long) :: fac2pi
  real(long) :: volg, volf
  real(long), pointer :: soug(:,:,:), souf(:,:,:)
!
  call alloc_array3d(soug, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(souf, 0, nrf, 0, ntf, 0, npf)
!
  omega = ome/radi
  call source_ang_mom_peos(souf)
  souf(0:nrf,0:ntf,0:npf) = omef(0:nrf,0:ntf,0:npf)*souf(0:nrf,0:ntf,0:npf)
  call vol_int_fluid(souf,volf)
  T_kinene_omeJ = 0.5d0*radi**3*volf
  W_gravene_omeJ = admmass - propermass - T_kinene_omeJ
!  T_kinene = 0.5d0*omega*angmom
!testtesst  call source_mp_minus_madm_peos(soug,souf)
!testtest  call vol_int_grav(soug,volg)
!testtest  call vol_int_fluid(souf,volf)
!
!testtest  fac2pi = 0.5d0/pi
!testtest  W_gravene = fac2pi*(radi*volg + radi**3*volf) + T_kinene
!
  ToverW_omeJ = T_kinene_omeJ/dabs(W_gravene_omeJ)
  I_inertia = angmom/omega
!
  deallocate(soug)
  deallocate(souf)
end subroutine calc_ToverW_peos
