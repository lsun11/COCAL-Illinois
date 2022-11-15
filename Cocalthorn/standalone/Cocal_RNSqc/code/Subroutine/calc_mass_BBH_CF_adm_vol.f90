subroutine calc_mass_BBH_CF_adm_vol
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg, rgin
  use make_array_2d
  use make_array_3d
  use def_quantities, only : admmass_thr
!  use interface_surf_source_adm_mass
  use interface_sourceterm_HaC_CF
  use interface_sourceterm_surface_int
  use interface_sourceterm_exsurf_eqm_binary
  use interface_bh_boundary_AH_n_psi
!  use coordinate_grav_r, only : rg
  use def_metric, only : psi
  use grid_parameter_binary_excision, only : ex_nrg
  use interface_vol_int_grav_bhex
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac2pi , fac4pi, val
  real(long) :: vol_int_adm, surf_int_bh, surf_int_ex
  real(long) :: adm_mass_vol, adm_mass_si_bh, adm_mass_si_ex
  real(long),pointer :: sou(:,:,:)   
  real(long), pointer :: sou_exsurf(:,:), dsou_exsurf(:,:)
  real(long), pointer :: sou_bhsurf(:,:), dsou_bhsurf(:,:)
  integer :: mass_ir, irg, ipg, itg
!
  call alloc_array3d(sou,0,nrg,0,ntg,0,npg)
  call alloc_array2d(sou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_bhsurf,0,ntg,0,npg)
  call alloc_array2d(sou_exsurf,0,ntg,0,npg)
  call alloc_array2d(dsou_exsurf,0,ntg,0,npg)

  fac2pi = 0.5d0/pi
  fac4pi = 0.25d0/pi

  call calc_mass_ir(mass_ir)
!  write(6,*) 'mass_ir = ', mass_ir

!  sou(0:nrg,0:ntg,0:npg) = 0.0d0
!  call excurve
  call sourceterm_HaC_CF(sou)
  call vol_int_grav_bhex(sou,vol_int_adm,mass_ir)
  adm_mass_vol = -fac2pi*vol_int_adm
!
  call sourceterm_surface_int(psi,0,sou_bhsurf,dsou_bhsurf)
  call surf_int_grav_rg(dsou_bhsurf,surf_int_bh,0)
  adm_mass_si_bh = -fac2pi*surf_int_bh
!
  call sourceterm_exsurf_eqm_binary(psi,sou_exsurf,dsou_exsurf)
  call surf_int_grav_rg(dsou_exsurf,surf_int_ex,ex_nrg)
  adm_mass_si_ex = -fac2pi*surf_int_ex
!
!  write(6,*) 'adm_mass_vol=', adm_mass_vol, 'adm_mass_si_bh=',adm_mass_si_bh, 'adm_mass_si_ex=',adm_mass_si_ex
  admmass_thr = adm_mass_vol + adm_mass_si_bh + adm_mass_si_ex
!
!  write (6,'(a14,1p,e14.6)') 'ADM mass thr =', admmass_thr
!
  dsou_bhsurf(0:ntg,0:npg) = 0.0d0
  call bh_boundary_AH_n_psi(dsou_bhsurf)
  call surf_int_grav_rg(dsou_bhsurf,surf_int_bh,0)
  adm_mass_si_bh = -fac2pi*surf_int_bh
!
  write (6,'(a40,1p,e20.12)') 'ADM mass thr using Neumann bc for psi=',  &
  &                          adm_mass_vol + adm_mass_si_bh + adm_mass_si_ex

  deallocate(sou)
  deallocate(sou_exsurf)
  deallocate(dsou_exsurf)
  deallocate(sou_bhsurf)
  deallocate(dsou_bhsurf)
end subroutine calc_mass_BBH_CF_adm_vol
