subroutine source_vep_surface_WL_peos(vpot_v,surp)
  use grid_parameter, only : nrf, ntf, npf
  use def_metric_on_SFC, only : hxxdf, hxydf, hxzdf, hyydf, hyzdf, hzzdf
  use def_matter, only : rs
  use def_velocity_potential, only : vep
  use interface_interpo_linear_surface_type0
  use interface_flgrad_midpoint_surface_type0
  implicit none
  real(long), pointer :: vpot_v(:,:,:), surp(:,:)
  real(long) :: hhxxu, hhxyu, hhxzu, hhyxu, hhyyu, hhyzu, &
  &             hhzxu, hhzyu, hhzzu
  real(long) :: gamxxu, gamxyu, gamxzu, gamyxu, gamyyu, gamyzu, &
  &             gamzxu, gamzyu, gamzzu
  real(long) :: dxvep, dyvep, dzvep
  integer    :: ir, it, ip
!
  ir = nrf
  do ip = 1, npf
    do it = 1, ntf
      call interpo_linear_surface_type0(hhxxu,hxxuf,ir,it,ip)
      call interpo_linear_surface_type0(hhxyu,hxyuf,ir,it,ip)
      call interpo_linear_surface_type0(hhxzu,hxzuf,ir,it,ip)
      call interpo_linear_surface_type0(hhyyu,hyyuf,ir,it,ip)
      call interpo_linear_surface_type0(hhyzu,hyzuf,ir,it,ip)
      call interpo_linear_surface_type0(hhzzu,hzzuf,ir,it,ip)
      hhyxu = hhxyu
      hhzxu = hhxzu
      hhzyu = hhyzu
      gamxxu = hhxxu + 1.0d0
      gamxyu = hhxyu
      gamxzu = hhxzu
      gamyyu = hhyyu + 1.0d0
      gamyzu = hhyzu
      gamzzu = hhzzu + 1.0d0
      gamyxu = gamxyu
      gamzxu = gamxzu
      gamzyu = gamyzu
!
      call flgrad_midpoint_surface_type0(vpot_v,dxvep,dyvep,dzvep,it,ip)
      call calc_surface_normal_midpoint(rs,rnx,rny,rnz,it,ip)
!
      surp(it,ip) = (gamxxu*dxvep + gamxyu*dyvep + gamxzu*dzvep)*rnx &
      &           + (gamyxu*dxvep + gamyyu*dyvep + gamyzu*dzvep)*rny &
      &           + (gamzxu*dxvep + gamzyu*dyvep + gamzzu*dzvep)*rnz
!
    end do
  end do
!
end subroutine source_vep_surface_WL_peos
