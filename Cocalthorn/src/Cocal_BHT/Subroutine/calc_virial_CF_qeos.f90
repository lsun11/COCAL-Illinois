subroutine calc_virial_CF_qeos
  use phys_constant, only  : long, pi
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use def_matter_parameter, only : radi
  use def_quantities, only : T_kinene, P_intene, W_gravene, &
  &                          ToverW, PoverW, Virial
  use make_array_3d
  use interface_source_virial_CF_qeos
  use interface_vol_int_grav
  use interface_vol_int_fluid
  implicit none
  real(long)          :: volf, volg, fac8pi
  real(long), pointer :: sou_Tkin(:,:,:), sou_Pint(:,:,:), sou_Wgra(:,:,:)
!
  call alloc_array3d(sou_Tkin, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(sou_Pint, 0, nrf, 0, ntf, 0, npf)
  call alloc_array3d(sou_Wgra, 0, nrg, 0, ntg, 0, npg)
!
  call source_virial_CF_qeos(sou_Tkin,sou_Pint,sou_Wgra)
!
  call vol_int_fluid(sou_Tkin,volf)
  T_kinene = radi**3*volf
!
  call vol_int_fluid(sou_Pint,volf)
  P_intene = radi**3*volf
!
  call vol_int_grav(sou_Wgra,volg)
  W_gravene= radi*volg
!
  ToverW = T_kinene/dabs(W_gravene)
  PoverW = P_intene/dabs(W_gravene)
  Virial = (2.0d0*T_kinene + 3.0d0*P_intene + W_gravene)/W_gravene
  Virial = dabs(Virial)
!  write (6,'(a20,1p,e14.6)') ' virial = ', virial
!
  deallocate(sou_Tkin)
  deallocate(sou_Pint)
  deallocate(sou_Wgra)
!
end subroutine calc_virial_CF_qeos
