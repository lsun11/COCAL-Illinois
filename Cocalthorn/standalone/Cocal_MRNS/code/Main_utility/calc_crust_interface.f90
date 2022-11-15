include '../Module/phys_constant.f90'
include '../EOS/Module/def_peos_parameter.f90'
include '../EOS/Subroutine/peos_initialize.f90'
include '../EOS/Subroutine/peos_lookup.f90'
!
program calc_crust_interface
!
  use phys_constant			!g,c,solmas,nnpeos
  use def_peos_parameter		!abc,abi,rhoi,qi,hi,nphase
  implicit none
!
  real(8) :: rho_int, abccgs_crust, abi_crust
!
  call peos_initialize
!
  abi_crust = 1.35692d+0
  abccgs_crust = 3.99874d-08*c**2
  rho_int = (abccgs(1)/abccgs_crust)**(1.0d0/(abi_crust - abi(1)))
  write(6,'(a23,es13.5)') 'crust Gamma            ', abi_crust
  write(6,'(a23,es13.5)') 'crust adiabtic constant', abccgs_crust
  write(6,'(a23,es13.5)') 'crust interface density', rho_int
!
end program calc_crust_interface
