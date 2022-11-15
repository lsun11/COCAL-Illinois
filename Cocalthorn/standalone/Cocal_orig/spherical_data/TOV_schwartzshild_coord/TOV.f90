include '../../code/Module/phys_constant.f90'
include '../../code/Module/def_matter_parameter.f90'
include '../../code/Module/def_quantities.f90'
include '../../code/Module/grid_parameter.f90'
include './Subroutine_TOV/rk.f90'
include './Subroutine_TOV/rkstep.f90'
include './Subroutine_TOV/equation_polytrope.f90'
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
  use def_quantities, only : restmass_sph, gravmass_sph, &
  &                          MoverR_sph, schwarz_radi_sph
  implicit none
!
  external equation
  integer, parameter :: neq = 6
  real(8)  :: y0(neq), yn(neq), work(neq, 2), &
           &  adm, compa, dqdq, dr, dr_back, &
           &  erer, ermx, radi_error, h, &
           &  compab, qini, qinib, radi, radib, radiini, &
           &  x0, xe, xn, yr, hini, ene, pre, rho0
  integer  :: i, ii, iter_compa, itmx, iter_radi, itype, idum, &
           &  ndiv, nstep
!
  ermx = 1.0d-9
  itmx = 1000
!
  open (2, file='ovpar.dat', status='old')
  open (1, file='ovphy.dat', status='unknown')
  open (8, file='ovlas.dat', status='unknown')
  open (9, file='ovaux.dat', status='unknown')
  read(2, '(2i5, 1es15.7)') nstep, ndiv, radiini
  read(2, '(2i5)') itype, idum
  close(2)
!
! --  read_parameter in grid_parameter.f90
!
  call read_parameter
!
! --  set compactness and initial
!
  qini = emdc
  dqdq = qini*0.1d0
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
      y0(1) = 0.0d+0 ; y0(2) = qini   ; y0(3) = 0.0d+0
      y0(4) = 0.0d+0 ; y0(5) = 2.0d+0 ; y0(6) = 0.0d+0
      yn(1) = 0.0d+0 ; yn(2) = qini   ; yn(3) = 0.0d+0 
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
        if (yn(2) <= 0.0d0) exit
!
!### need to be changed for output profiles
!        if (radi_error <= ermx .and. itype == 1) then
!          write(8, '(7(es15.7))') xe, yn(2), yn(2)**pinx, yn(2)**(1.0+pinx), yn(5)
!        end if
      end do
!
      if (radi_error <= ermx) then
        if (itype == 0) exit
        if (itype == 1) then
          yr = xe/yn(1)*(1.0d0-(1.0d0-2.0d0*yn(1)/xe)**0.5)
          adm =  yn(5)*yn(6)/yr
          compa = yn(1)/xe
          write(6, '(5es14.6)') xe, yn(1), yn(2), yn(3), yn(4), &
          &                     yn(5), yn(6), qini, compa, adm
          write(1, '(5es14.6)') xe, yn(1), yn(2), yn(3), yn(4), &
          &                     yn(5), yn(6), qini, compa, adm
          stop ' end execution '
        end if
      end if
!
! ---- hunting for precise radius.
!
      if (yn(2) <= 0.0d0) then
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
    if (iter_compa == 1) then
      if (compa > MoverR_sph) dqdq = - dqdq
      compab = compa
      qinib= qini
      qini = qini + dqdq
      dqdq = dabs(dqdq)
      erer = 1.0d0
      write(6, '(1i5, es12.4, 3es16.8)') iter_compa, erer, qinib, compa
      cycle
    end if
!
    dqdq = - (compa-MoverR_sph)/(compa-compab)*(qini-qinib)
    compab = compa
    qinib= qini
    qini = qini + dqdq
!
    erer = dabs((compa-MoverR_sph)/MoverR_sph)
    write(6, '(1i5, es12.4, 3es16.8)') iter_compa, erer, qinib, compa
    if (erer < ermx) then
      write(6 , '(5es14.6)') xe, yn(1), yn(2), yn(3), yn(4), &
      &                      yn(5),yn(6),qinib, compa, adm
      write(1 , '(5es14.6)') xe, yn(1), yn(2), yn(3), yn(4), &
      &                      yn(5),yn(6),qinib, compa, adm
!
      open(14, file='rnspar_add.dat', status='unknown') 
!
      emdc = qini
      restmass_sph = yn(3)
      gravmass_sph = yn(1)
      MoverR_sph  = compa
      schwarz_radi_sph = xe
      write(14, '(2es14.6, a15)') emdc, pinx, &
        & '  : emdc, pinx'
      write(14, '(2es14.6, a31)') restmass_sph, gravmass_sph, &
        & '  : restmass_sph, gravmass_sph'
      write(14, '(2es14.6, a47)') MoverR_sph, schwarz_radi_sph, &
        & '  : MoverR_sph,  schwarz_radi_sph  (K=1 unit)'
      write(6, '(/, /, 2es14.6, a15)') emdc, pinx, &
        & '  : emdc, pinx'
      write(6, '(2es14.6, a31)') restmass_sph, gravmass_sph, &
        & '  : restmass_sph, gravmass_sph'
      write(6, '(2es14.6, a47)') MoverR_sph, schwarz_radi_sph, &
        & '  : MoverR_sph,   schwarz_radi_sph  (K=1 unit)'
!
      close(14)
!
      stop ' end of execution '
    end if
!
  end do
!
  stop
end PROGRAM TOV
