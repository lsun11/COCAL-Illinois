! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
! G = c = Msol = 1 unit is used in TOV integration.
! Use pre_c dyn/cm^2 at rho_c gr/cm^3 to specify K.
! pre_0 and rho_0 are input parameter,  
! typically pre_0 = 1.0d+37 dyn/cm^2 
! and       rho_0 = 1.0d+16  gr/cm^3.
!
! --- Desctiption for input parameters.
!
!  -- Input in 'ovpar_parameos.dat'
!
! hini  : not used.
! nstep : see Runge Kutta subroutine.
! ndiv  : see Runge Kutta subroutine.
! radini: initial radius for the surface of a neutron star.
!        (Iteration is made until "radi" to become a surface exactly.)
! itype : itype = 0  iteration for compactness
!         itype = 1  single structure
! chope : compactness to find (for itype = 0)
!
!
!  -- Input in "oveosparm.dat"
!
! nphase : total number of piecewise regions.
! rhoini : central value of the density [g/cm^3]. 
!         (initial for the integration.)
! rhocgs : value of the density in a piecewise region.  
! abi    : value of the adiabatic index (Gamma) in a piecewise region.  
!
!  -- Ouput in "oveos_phase.dat" : data at the interface of each piecewise 
!                              region, rescaled in G=c=Msol=1 unit.
!
!  write(30,*)'#', nphase, qini
!  write(30,40) ii, abc(ii), abi(ii), rhoi(ii), qi(ii), hi(ii), 
! &              rhocgs(ii), abccgs(ii)
!
!  qini   : p/rho at the center rescaled in G=c=Msol=1 unit.
!  abc    : adiabatic constant K_i rescaled in G=c=Msol=1 unit.
!  abi    : adiabatic index Gamma_i.
!  rhoi   : rest mass density rescaled in G=c=Msol=1 unit.
!  qi     : p/rho rescaled in G=c=Msol=1 unit.
!  hi     : the relativistic enthalpy rescaled in G=c=Msol=1 unit.
!  rhocgs : rest mass density in [g/cm^3].
!  abccgs : adiabatic constant K_i rescaled in [g/cm^3].
! ______________________________________________________________________
! ______________________________________________________________________
!
include './Module/phys_constant.f90'
include './Module/def_peos_parameter.f90'
include './Subroutine/peos_initialize.f90'
include './Subroutine/peos_lookup.f90'
include './Subroutine/peos_q2hprho.f90'
include './Subroutine/peos_h2qprho.f90'
!
PROGRAM main
  implicit none
  real(8) :: q, h, pre, rho, ened
!
  call peos_initialize
!
  q = 2.25576d-01
  call peos_q2hprho(q,h,pre,rho)
  open(870,file='peos_par.dat',status='unknown')
  write(870,'(4es13.5, a28)') q, h, pre, rho, ' : q, h, pre, rho'
!
  call peos_h2qprho(h,q,pre,rho,ened)
  write(870,'(5es13.5, a34)') h, q, pre, rho, ened, ' : h, q, pre, rho, ened'
  close(870)
!
end PROGRAM main
