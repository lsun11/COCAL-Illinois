include '../../code/Module/phys_constant.f90'
include '../../code/Module/def_matter.f90'
include '../../code/Module/def_matter_parameter.f90'
include '../../code/Module/def_quantities.f90'
include '../../code/Module/grid_parameter.f90'
include '../../code/EOS/Module/def_peos_parameter.f90'
include '../../code/EOS/Subroutine/peos_initialize.f90'
include '../../code/EOS/Subroutine/peos_lookup.f90'
include '../../code/EOS/Subroutine/peos_h2qprho.f90'
include '../../code/EOS/Subroutine/peos_q2hprho.f90'
include './Module_TOV/def_TOV_quantities.f90'
include './Subroutine_TOV/rk.f90'
include './Subroutine_TOV/rkstep.f90'
include './Subroutine_TOV/equation_peos.f90'
include './Subroutine_TOV/printout_TOV_profile.f90'
include './Subroutine_TOV/printout_TOV_quantities.f90'
!
PROGRAM TOV
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
! --- Desctiption for parameters.
!
!    qini  : initial for emden function. 
!    pinx  : polytropic index
!    nstep : see Runge Kutta subroutine.
!    ndiv  : see Runge Kutta subroutine.
!    radini: initial radius for hunting radius from inside.
!    itype : itype = 0  iteration for compactness
!            itype = 1  single structure
!    chope : compactness to find (for itype = 0)
!
! --- Description for variables.
!
!    xe   : schwartzscild radial coordinate 
!    yn(1) : effective gravitational mass
!    yn(2) : emden function
!    yn(3) : proper rest mass 
!    yn(4) : proper mass energy
!    yn(5) : conformal factor (proper boundary condition not applied)
!    yn(6) : ADM mass energy  (proper boundary condition not applied)
!    adm  : ADM mass energy
!    compa : Compaction = YN(1)/XE = 
!              = Effective gravitational mass/
!                Radius of star in Schwartzscild radial coordinate 
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
  use grid_parameter
  use def_matter_parameter, only : emdc, pinx
  use def_peos_parameter, only : emdini_gcm1
  use def_quantities, only : restmass_sph, gravmass_sph, &
  &                          MoverR_sph, schwarz_radi_sph
  use def_TOV_quantities
  use phys_constant, only: g, c, solmas
  implicit none
!
  external equation
  real(8)  :: Q_target, Q_current, Q_current_bak
  integer  :: i, ii, iter_compa, itmx, iter_radi, itype, iter_mode, &
           &  ndiv, nstep
  character(2) :: file_ctl
!
  ermx      = 1.0d-10
  ermx_MorC = 1.0d-9
  itmx = 1000
!
  open (2, file='ovpar_peos.dat', status='old')
  open (9, file='ovaux.dat', status='unknown')
  read(2, '(2i5, 1es15.7)') nstep, ndiv, radiini
  read(2, '(2i5)') itype, iter_mode
  close(2)
!
! --  read_parameter in grid_parameter.f90
! --  set compactness and initial
!
  call read_parameter
!
  call peos_initialize
  qini = emdini_gcm1
  call peos_q2hprho(qini, hini, pre, rho0, ene)
!  dqdq = hini*0.001d0
  dqdq = hini*0.01d0
!  dqdq = hini*0.1d0
!
! --- Iteration for compactness.
!
  do iter_compa = 1, itmx
!
! --- Iteration for radius.
!
    radi = radiini
    dr = 0.2d0*radi
    dr_back = radi
    radi_error = dr/radi
!
    do iter_radi = 1, itmx
!
      x0 = 0.0d+0 ; xn = radi
      y0(1) = 0.0d+0 ; y0(2) = hini   ; y0(3) = 0.0d+0
      y0(4) = 0.0d+0 ; y0(5) = 2.0d+0 ; y0(6) = 0.0d+0
      yn(1) = 0.0d+0 ; yn(2) = hini   ; yn(3) = 0.0d+0 
      yn(4) = 0.0d+0 ; yn(5) = 2.0d+0 ; yn(6) = 0.0d+0
      h = (xn - x0) / dble(nstep)
!
!### need to be changed for output profiles
!      if (radi_error <= ermx .and. itype == 1) &
!      & write(8 , '(7(es15.7))') xe, yn(2), yn(2)**pinx, yn(2)**(1.0d0+pinx), yn(5)
!
      do i = 1, nstep
        ii = i
        xe = x0 + h
        call rk(neq, equation, x0, xe, ndiv, y0, yn, work)
        if (yn(2) <= 1.0d0) exit
!
!### need to be changed for output profiles
        if (radi_error <= ermx .and. itype == 1) then
          file_ctl = 'nm'
          if (i.eq.1)     file_ctl = 'op'
          if (i.eq.nstep) file_ctl = 'cl'
          call printout_TOV_profile(file_ctl)
        end if
      end do
!
      if (radi_error <= ermx) then
        if (itype == 0) exit
        if (itype == 1) then
          call printout_TOV_quantities
          stop
!          stop ' end execution '
        end if
      end if
!
! ---- hunting for precise radius.
!
      if (yn(2) <= 1.0d0) then
        if (iter_radi == 1) write(6, *) ' bad initial radius '
        radi = radi - dr
        dr_back = dr
        write(9, *)' back ', iter_radi, ii, radi, radi_error
      else
        radi = radi + 0.2*dr_back
        radi_error = dr_back/radi
        dr = 0.2d0*dr_back
        write(9, *)' hunt ', iter_radi, ii, radi, radi_error
      end if
      if (iter_radi == itmx) stop ' iter_radi '
!
! ---- hunting end.
!
    end do
!
! ---- Searching for a certain value of compactness.
!
    yr = xe/yn(1)*(1.0d0-(1.0d0-2.0d0*yn(1)/xe)**0.5)
    adm = yn(5)*yn(6)/yr
    compa = yn(1)/xe
!
    if (iter_mode.eq.0) then
      Q_current = compa
      Q_target  = MoverR_sph
    else if (iter_mode.eq.1) then
      Q_current = yn(3)
      Q_target  = restmass_sph
    else if (iter_mode.eq.2) then
      Q_current = yn(1)
      Q_target  = gravmass_sph
    else 
      stop 'set iteration mode for a quantity'
    end if
!
    if (iter_compa == 1) then
      if (Q_current > Q_target) dqdq = - dqdq
      Q_current_bak = Q_current
      hinib= hini
      hini = hini + dqdq
      dqdq = dabs(dqdq)
      erer = 1.0d0
      write(6, '(1i5, es12.4, 3es16.8)') iter_compa, erer, hinib, Q_current
      cycle
    end if
!
    dqdq = - (Q_current-Q_target)/(Q_current-Q_current_bak)*(hini-hinib)
    Q_current_bak = Q_current
    hinib= hini
    hini = hini + dqdq
!
    erer = dabs((Q_current-Q_target)/Q_target)
    write(6, '(1i5, es12.4, 3es16.8)') iter_compa, erer, hinib, Q_current
!
!
    if (erer < ermx_MorC) then
      call printout_TOV_quantities
!
!      stop ' end of execution '
      stop
    end if
!
  end do
!
  stop
end PROGRAM TOV
