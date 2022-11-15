subroutine interpolation_cartesian_MRNS
  use phys_constant, only : long
  use grid_parameter, only : nrg, ntg, npg, nrf, ntf, npf
  use grid_parameter_cartesian, only : nx, ny, nz
  use coordinate_grav_r, only : rg
  use trigonometry_grav_theta, only : sinthg
  use trigonometry_grav_phi, only : sinphig, cosphig
  use def_metric
  use def_metric_hij
  use def_emfield
  use def_faraday_tensor
  use def_matter, only : emd, utf, uxf, uyf, uzf
  use def_metric_cartesian
  use def_metric_hij_cartesian
  use def_emfield_cartesian
  use def_faraday_tensor_cartesian
  use def_matter_cartesian
  use def_matter_velocity
  use interface_modules_cartesian
  use make_array_3d
  implicit none
  real(long) :: xxxx, yyyy
  real(long), pointer :: fnc_in(:,:,:), fnc_out(:,:,:)
  integer :: ir, it, ip
!
  call alloc_array3d(fnc_in,0,nrg,0,ntg,0,npg)
  call alloc_array3d(fnc_out,1,nx,1,ny,1,nz)
!
  call interpolation_metric(psi,psica)
  call interpolation_metric(alph,alphca)
  call interpolation_metric(bvxd,bvxdca)
  call interpolation_metric(bvyd,bvydca)
  call interpolation_metric(bvzd,bvzdca)
  call interpolation_metric(hxxd,hxxdca)
  call interpolation_metric(hxyd,hxydca)
  call interpolation_metric(hxzd,hxzdca)
  call interpolation_metric(hyyd,hyydca)
  call interpolation_metric(hyzd,hyzdca)
  call interpolation_metric(hzzd,hzzdca)
!
call invhij
call calc_shift_down2up
  call interpolation_metric(bvxu,bvxuca)
  call interpolation_metric(bvyu,bvyuca)
  call interpolation_metric(bvzu,bvzuca)
!
  call interpolation_metric(va,vaca)
  call interpolation_metric(vaxd,vaxdca)
  call interpolation_metric(vayd,vaydca)
  call interpolation_metric(vazd,vazdca)
  call interpolation_metric(fxd_grid,fxd_gridca)
  call interpolation_metric(fyd_grid,fyd_gridca)
  call interpolation_metric(fzd_grid,fzd_gridca)
  fnc_in(0:nrg,0:ntg,0:npg)     = fijd_grid(0:nrg,0:ntg,0:npg,1)
  call interpolation_metric(fnc_in,fnc_out)
  fijd_gridca(1:nx,1:ny,1:nz,1) = fnc_out(1:nx,1:ny,1:nz)
  fnc_in(0:nrg,0:ntg,0:npg)     = fijd_grid(0:nrg,0:ntg,0:npg,2)
  call interpolation_metric(fnc_in,fnc_out)
  fijd_gridca(1:nx,1:ny,1:nz,2) = fnc_out(1:nx,1:ny,1:nz)
  fnc_in(0:nrg,0:ntg,0:npg)     = fijd_grid(0:nrg,0:ntg,0:npg,3)
  call interpolation_metric(fnc_in,fnc_out)
  fijd_gridca(1:nx,1:ny,1:nz,3) = fnc_out(1:nx,1:ny,1:nz)
!
  call interpolation_matter(emd,emdca)
  call interpolation_matter(utf,utca)
  vxu(0:nrf,0:ntf,0:npf) = uxf(0:nrf,0:ntf,0:npf)/utf(0:nrf,0:ntf,0:npf)
  vyu(0:nrf,0:ntf,0:npf) = uyf(0:nrf,0:ntf,0:npf)/utf(0:nrf,0:ntf,0:npf)
  vzu(0:nrf,0:ntf,0:npf) = uzf(0:nrf,0:ntf,0:npf)/utf(0:nrf,0:ntf,0:npf)
  call interpolation_matter(vxu,vxca)
  call interpolation_matter(vyu,vyca)
  call interpolation_matter(vzu,vzca)
!
  deallocate(fnc_in)
  deallocate(fnc_out)
end subroutine interpolation_cartesian_MRNS
