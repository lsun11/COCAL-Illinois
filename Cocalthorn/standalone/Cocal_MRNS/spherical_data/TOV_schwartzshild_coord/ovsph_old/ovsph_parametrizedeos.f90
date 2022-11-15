include '../../code/EOS/Module/phys_constant.f90'
include '../../code/Module/def_quantities.f90'
include '../../code/Module/def_matter_parameter.f90'
include '../../code/Module/grid_parameter.f90'
include '../../code/EOS/Module/def_peos_parameter.f90'
include '../../code/EOS/Subroutine/peos_initialize.f90'
include '../../code/EOS/Subroutine/peos_lookup.f90'
include '../../code/EOS/Subroutine/peos_q2hprho.f90'
include '../../code/EOS/Subroutine/peos_h2qprho.f90'
include './Subroutine_TOV/rk.f90'
include './Subroutine_TOV/rkstep.f90'
include './Subroutine_TOV/equation_2.f90'
!
PROGRAM TOV
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
!     G = c = Msol = 1 unit is used in TOV integration.
!     assume pre = 10^37 dyn/cm^2 at rho = 10^16 gr/cm^3 to specify K.
!
! --- Desctiption for input parameters.
!
!  -- Input in 'ovpar_parameos.dat'
!
!     hini  : not used.
!     nstep : see Runge Kutta subroutine.
!     ndiv  : see Runge Kutta subroutine.
!     radini: initial radius for the surface of a neutron star.
!            (Iteration is made until "radi" to become a surface exactly.)
!     itype : itype = 0  iteration for compactness
!             itype = 1  single structure
!     chope : compactness to find (for itype = 0)
!
!
!  -- Input in "oveosparm.dat"
!
!     nphase : total number of piecewise regions.
!     rhoini : central value of the density [g/cm^3]. 
!             (initial for the integration.)
!     rhocgs : value of the density in a piecewise region.  
!     abi    : value of the adiabatic index (Gamma) in a piecewise region.  
!
! --- Description for the variables.
!
!  -- Output in "ovphy.dat" : Mass, radius and central enthalpy
!
!       write(1, 630) xe, yn(1), yn(2), yn(3), yn(4), 
!     &                 yn(5), yn(6), hini, compa, adm
!
!     xe    : schwartzscild radial coordinate at the surface
!     yn(1) : gravitational mass
!     yn(2) : pressure
!     yn(3) : proper rest mass 
!     yn(4) : proper mass energy
!     yn(5) : conformal factor (proper boundary condition is not applied)
!     yn(6) : qantity relates to ADM mass energy.  
!     hini  : the relativistic enthalpy at the center.
!     adm   : ADM mass energy
!     compa : Compactness = YN(1)/XE = adm/XE
!               = Gravitational mass/
!                 Radius of star in Schwartzscild radial coordinate 
!
!  -- Output in "ovlas.dat" : profile of thermodynamic variables.
!
!     write(8, 600) xe, hh, rho0, pre, ene, yn(5)
!
!     xe    : schwartzscild radial coordinate 
!     hh    : the relativistic enthalpy.
!     rho0  : the rest mass density.
!     pre   : the pressure.
!     ene   : the energy density.
!     yn(5) : conformal factor (proper boundary condition is NOT applied)
!
!  -- Ouput in "oveos_phase.dat" : data at the interface of each piecewise 
!                                  region, rescaled in G=c=Msol=1 unit.
!
!      write(30, *)'#', nphase, qini
!      write(30, 40) ii, abc(ii), abi(ii), rhoi(ii), qi(ii), hi(ii), 
!     &              rhocgs(ii)
!
!      qini   : p/rho at the center rescaled in G=c=Msol=1 unit.
!      abc    : adiabatic constant K_i rescaled in G=c=Msol=1 unit.
!      abi    : adiabatic index Gamma_i.
!      rhoi   : rest mass density rescaled in G=c=Msol=1 unit.
!      qi     : p/rho rescaled in G=c=Msol=1 unit.
!      hi     : the relativistic enthalpy rescaled in G=c=Msol=1 unit.
!      rhocgs : rest mass density in [g/cm^3].
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
  use phys_constant, only : pi
  use grid_parameter, only : pinx
  use def_peos_parameter, only : abc, abi, qi, hi, nphase
  implicit none
!
  external equation
  integer, parameter :: neq = 6
  real(8) :: y0(neq), yn(neq), work(neq, 2)
  real(8) :: abct, abin, adm, chope, compa, compab, dqdq, dr, &
  &          ene, erer, ermx, radi_error, h, hh, hini, hinib, hinibx, &
  &          pre, q, qini, radi, radib, radiini, rho0, x0, &
  &          xe, xn, yr
  integer :: i, idum, ii, iph, iter_compa, itmx, iter_radi, itype, ndiv, &
  &          nstep
!
  ermx = 1.0d-10
  itmx = 100
  qini = 5.81922018
!
  open (2, file='ovpar_parameos.dat', status='old')
  open (1, file='ovphy.dat', status='unknown')
  open (8, file='ovlas.dat', status='unknown')
  open (9, file='ovaux.dat', status='unknown')
  read(2, '(2i5, 1es15.7)') nstep, ndiv, radiini
  read(2, '(2i5, 1es15.7)') itype, idum, chope
  close(2)
!
  call peos_initialize
  call peos_q2hprho(qini, hini, pre, rho0, ene)
!
!  hh = hini
!
! --- iteration for compactness.
!
  dqdq = hini*0.1d0
!
  do iter_compa = 1, itmx
!
! --- iteration for radius.
!
    radi = radiini
    dr = 0.2d0*radi
    radi_error = dr/radi
!
    do iter_radi = 1, itmx
!
      x0 = 0.0d+0 ; xn = radi
      y0(1) = 0.0d+0 ; y0(2) = hini ; y0(3) = 0.0d+0 ; y0(4) = 0.0d+0
      y0(5) = 2.0d+0 ; y0(6) = 0.0d+0
      yn(1) = 0.0d+0 ; yn(2) = hini ; yn(3) = 0.0d+0 ; yn(4) = 0.0d+0
      yn(5) = 2.0d+0 ; yn(6) = 0.0d+0
      h = (xn - x0) / nstep
!
      if (radi_error <= ermx.and.itype == 1) then
        call peos_h2qprho(hini, q, pre, rho0, ene)
        write(8, '(7(es15.7))') x0, hini, rho0, pre, ene, yn(5)
      end if
!
      do i = 1, nstep
        ii = i
        xe = x0 + h
        call rk(neq, equation, x0, xe, ndiv, y0, yn, work)
!
        hh = yn(2)
        call peos_h2qprho(yn(2), q, pre, rho0, ene)
!
        if (radi_error <= ermx .and. itype == 1) &
        & write(8, '(7(es15.7))') xe, yn(2), rho0, pre, ene, yn(5)
!
        if (yn(2) <= 1.0d0) exit
!
      end do
!
      if (radi_error <= ermx) then
        if (itype == 0) exit
        if (itype == 1) then
          yr = xe/yn(1)*(1.0d0-(1.0d0-2.0d0*yn(1)/xe)**0.5)
          adm =  yn(5)*yn(6)/yr
          compa = yn(1)/xe
          write(6, '(5es14.5)') xe, yn(1), yn(2), yn(3), yn(4), &
            &                 yn(5), yn(6), hini, compa, adm
          write(1, '(5es14.5)') xe, yn(1), yn(2), yn(3), yn(4), &
            &                 yn(5), yn(6), hini, compa, adm
          stop ' end execution '
        end if
      end if
!
! ---- look for a precise radius.
!
! --  using hh
!
      if (yn(2) <= 1.0d0) then
        if (iter_radi == 1) then
          write(6, *) ' bad initial radius '
          radi = 0.95d0*radi
        else
          radi = radib
          radi_error = dr/radi
          dr = 0.2d0*dr
          write(9, *)' back ', iter_radi, ii, radi, radi_error
        end if
      else
        radib = radi
        radi = radi + dr
        write(9, *)' hunt ', iter_radi, ii, radi, radi_error, yn(2)
      end if
!
      if (iter_radi == itmx) stop ' iter '
    end do
!
! ---- look for the radius end.
!
! ---- searching for a certain value of compactness.
!
    yr = xe/yn(1)*(1.0d0-(1.0d0-2.0d0*yn(1)/xe)**0.5)
    adm =  yn(5)*yn(6)/yr
    compa = yn(1)/xe
!
    if (iter_compa == 1) then
      if (compa > chope) dqdq = - dqdq
      compab = compa
      hinib= hini
      hini = hini + dqdq
      dqdq = dabs(dqdq)
      erer = 1.0d0
      write(6, '(1i5, es12.4, 3es16.8)') iter_compa, erer, hinib, compa
      cycle
    end if
!
    if (chope >= compa.and.compa >= compab) then
      compab = compa
      hinib= hini
      hini = hini + dqdq
    else if (chope <= compa.and.compa <= compab) then
      compab = compa
      hinib= hini
      hini = hini - dqdq
    else
      compab = compa
      hinibx = hini
      hini = 0.5d0*(hini+hinib)
      hinib = hinibx
      dqdq = 0.5d0*dabs(dqdq)
    end if
!
    erer = dabs((compa-chope)/chope)
    write(6, '(1i5, es12.4, 3es16.8)') iter_compa, erer, hinib, compa
    if (erer < ermx) then
      write(6 , '(5es14.5)') xe, yn(1), yn(2), yn(3), yn(4), &
        &                    yn(5), yn(6), hinib, compa, adm
      write(1 , '(5es14.5)') xe, yn(1), yn(2), yn(3), yn(4), &
        &                    yn(5), yn(6), hinib, compa, adm
      stop ' end of execution '
    end if
!
  end do
!
  stop 
END PROGRAM TOV
