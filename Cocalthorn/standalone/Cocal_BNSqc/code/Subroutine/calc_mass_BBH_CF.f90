subroutine calc_mass_BBH_CF
  use phys_constant, only  :   long, pi
  use grid_parameter, only :   nrg, ntg, npg
  use make_array_1d
  use make_array_2d
  use def_quantities, only : admmass, komarmass, admmass_asymp, komarmass_asymp
  use interface_surf_source_adm_mass
  use interface_surf_source_komar_mass
!  use interface_source_komar_mass_compact
!  use interface_surf_int_grav
  use interface_surf_int_grav_rg
  implicit none
  real(long) :: fac2pi , fac4pi
  real(long) :: surf_int_adm, surf_int_kom
  real(long) :: adm_em, komar_em
  real(long),pointer :: surf_sou_adm(:,:), surf_sou_kom(:,:)
  real(long),pointer :: adm_mass(:), komar_mass(:)
  integer :: mass_ir, ii, irg, Nmass
!
  Nmass = 8
  call alloc_array1d(adm_mass  ,1,Nmass)
  call alloc_array1d(komar_mass,1,Nmass)
  adm_mass(1:Nmass)  =0.0d0
  komar_mass(1:Nmass)=0.0d0
!
  call alloc_array2d(surf_sou_adm, 1, ntg, 1, npg)
  call alloc_array2d(surf_sou_kom, 1, ntg, 1, npg)
!
  fac2pi = 0.5d0/pi
  fac4pi = 0.25d0/pi
  call calc_mass_ir(mass_ir)
  write(6,'(a10,i5)') 'mass_ir = ', mass_ir
!
  do ii=1,Nmass
    surf_sou_adm(1:ntg,1:npg) = 0.0d0
    surf_sou_kom(1:ntg,1:npg) = 0.0d0
    surf_int_adm = 0.0d0
    surf_int_kom = 0.0d0
    adm_em = 0.0d0
    komar_em = 0.0d0
!
    irg = mass_ir+ii-Nmass/2
!
    call surf_source_adm_mass(surf_sou_adm,irg)
    call surf_int_grav_rg(surf_sou_adm, surf_int_adm, irg)
    adm_mass(ii) = -fac2pi*surf_int_adm
!
    call surf_source_komar_mass(surf_sou_kom, irg)
    call surf_int_grav_rg(surf_sou_kom, surf_int_kom, irg) 
    komar_mass(ii) = fac4pi*surf_int_kom
!    write (6,'(a6,i5,a20,1p,e14.6)') 'irg = ', irg, '   ADM   mass =     ', adm_mass(ii)
!    write (6,'(a6,i5,a20,1p,e14.6)') 'irg = ', irg, '   Komar mass =     ', komar_mass(ii)
!
    if (ii.ge.2) then
      adm_em = dabs(adm_mass(ii)-adm_mass(ii-1))/adm_mass(ii)
      komar_em = dabs(komar_mass(ii)-komar_mass(ii-1))/komar_mass(ii)
    end if
    write(6,'(a6,i5,a15,1p,e14.6,a15,1p,e14.6,a15,1p,e14.6,a15,1p,e14.6)') 'irg = ',irg,   &
      & '    ADM mass = ',adm_mass(ii), '  Komar mass = ',komar_mass(ii), &
      & '    ADM error= ',adm_em,       '  Komar error= ',komar_em
  end do  
  admmass = adm_mass(Nmass/2)
  komarmass = komar_mass(Nmass/2)

  admmass_asymp = admmass
  komarmass_asymp = komarmass

  deallocate(adm_mass)
  deallocate(komar_mass)
  deallocate(surf_sou_adm)
  deallocate(surf_sou_kom)
end subroutine calc_mass_BBH_CF
