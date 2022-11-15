subroutine calc_EMenergy_axisym_WL
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_matter_parameter, only : radi
  use def_quantities, only : M_torBene, M_polBene, M_eleEene, W_gravene, &
  &                          MtorBoverW, MpolBoverW, MeleEoverW
  use make_array_3d
  use interface_source_EMenergy_axisym_WL
  use interface_vol_int_grav
  implicit none
  real(long)          :: volg, fac8pi
  real(long), pointer :: sou_MtorB(:,:,:), sou_MpolB(:,:,:), &
  &                      sou_MeleE(:,:,:)
!
  call alloc_array3d(sou_MtorB, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(sou_MpolB, 0, nrg, 0, ntg, 0, npg)
  call alloc_array3d(sou_MeleE, 0, nrg, 0, ntg, 0, npg)
!
  call source_EMenergy_axisym_WL(sou_MtorB,sou_MpolB,sou_MeleE)
!
  call vol_int_grav(sou_MtorB,volg)
  M_torBene = radi*volg
!
  call vol_int_grav(sou_MpolB,volg)
  M_polBene = radi*volg
!
  call vol_int_grav(sou_MeleE,volg)
  M_eleEene = radi*volg
!
  MtorBoverW = M_torBene/dabs(W_gravene)
  MpolBoverW = M_polBene/dabs(W_gravene)
  MeleEoverW = M_eleEene/dabs(W_gravene)
!
  deallocate(sou_MtorB)
  deallocate(sou_MpolB)
  deallocate(sou_MeleE)
!
end subroutine calc_EMenergy_axisym_WL
