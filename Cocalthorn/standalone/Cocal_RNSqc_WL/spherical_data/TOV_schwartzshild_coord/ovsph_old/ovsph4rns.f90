include './Module/phys_constant.f90'
include './Module/def_matter.f90'
include './Module/def_quantities.f90'
include './Module/grid_parameter.f90'
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
  use grid_parameter	!pinx, emdc, restmass_sph, gravmass_sph, MoverR_sph, schwarz_radi_sph
  implicit none
!
  external test
  integer, parameter :: neq = 6
  real(8)  :: y0(neq), yn(neq), work(neq, 2), &
           &  adm, compa, dqdq, dr, dum, &
           &  erer, ermx, radi_error, h, &
           &  compab, pini, qini, qinib, radi, radib, radiini, &
           &  x0, xe, xn, yr
  integer  :: i, ii, iter_compa, itmx, iter_radi, itype, idum, &
           &  ndiv, ndum, nstep
!
  ermx = 1.0d-9
  itmx = 100
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
  call read_parameter	!emdc, restmass_sph, gravmass_sph, MoverR_sph, schwarz_radi_sph
!
! --  set compactness and initial
!
  qini = emdc
!
! --- Iteration for compactness.
!
  dqdq = qini*0.1d0
!
  do iter_compa = 1, itmx	!2222
!
! --- Iteration for radius.
!
    radi = radiini
    dr = 0.2d0*radi
    radi_error = dr/radi
!
    do iter_radi = 1, itmx	!200
!
      pini = qini**(pinx+1.0d0)
      x0 = 0.0d+0 ; xn = radi
      y0(1) = 0.0d+0 ; y0(2) = qini ; y0(3) = 0.0d+0 ; y0(4) = 0.0d+0 ; y0(5) = 2.0d+0 ; y0(6) = 0.0d+0
      yn(1) = 0.0d+0 ; yn(2) = qini ; yn(3) = 0.0d+0 ; yn(4) = 0.0d+0 ; yn(5) = 2.0d+0 ; yn(6) = 0.0d+0
      h = (xn - x0) / nstep
!
      if (radi_error <= ermx .and. itype == 1) &
      & write(8 , '(7(es15.7))') xe, yn(2), yn(2)**pinx, yn(2)**(1.0+pinx), yn(5)
!
      do i = 1, nstep
        ii = i
        xe = x0 + h
        call rk(neq, test, x0, xe, ndiv, y0, yn, work)
        if (yn(2) <= 0.0d0) exit
        if (radi_error <= ermx .and. itype == 1) &
        & write(8, '(7(es15.7))') xe, yn(2), yn(2)**pinx, yn(2)**(1.0+pinx), yn(5)
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
        if (iter_radi == 1) then
          write(6, *) ' bad initial radius '
          radib = 0.90d0*radi					!™
          radi = 0.95d0*radi
        else
          radi = radib
          radi_error = dr/radi
          dr = 0.2d0*dr
!          radi_error = dr/radib			!š
!          dr = 0.2d0*dr				!š
!          radi = radib + dr				!š
          write(9, *)' back ', iter_radi, ii, radi, radi_error
        end if
      else
        radib = radi
        radi = radi + dr
        write(9, *)' hunt ', iter_radi, ii, radi, radi_error
      end if
!
! ---- hunting end.
!
    end do	!200
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
      cycle	!2222
    end if
!
    dqdq = - (compa-MoverR_sph)/(compa-compab)*(qini-qinib)
    compab = compa
    qinib= qini
    qini = qini + dqdq
!
    erer = dabs((compa-MoverR_sph)/MoverR_sph)
!  write(6, *) compab, compa, dqdq
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
      moverr_sph  = compa
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
  end do	!2222
!
  stop
end program tov
!
! --------------------------------------------------------
! --- equations are in this routine.
!
subroutine test(x, y, f)
  use def_matter, only : pinx
  use phys_constant, only : pi
  implicit none
!
  integer, parameter     :: neq = 6
  real(8), intent(inout) :: y(neq), f(neq), x
  real(8) :: rr, rma, emd, bma, tma, psi, adm, &
          &  rho0, rho, pre, hhh
!
  rr  = x
  rma = y(1)
  emd = y(2)
  bma = y(3)
  tma = y(4)
  psi = y(5)
  adm = y(6)
  if (emd <= 0.0d0) then
    rho0= 0.0d0
    rho = 0.0d0
    pre = 0.0d0
    hhh = 1.0d0
  else
    rho0= emd**pinx
    rho = emd**pinx*(1.0d0 + pinx*emd)
    pre = emd**(1.0d0 + pinx)
    hhh = 1.0d0 + (1.0d0+pinx)*emd
  end if
!
  if (rr == 0.0d0) then
    f(1:6) = 0.0d0
  else
    f(1) = 4.0d0*pi*rr**2*rho
    f(2) = - 1.0d0/(1.0d0+pinx)*hhh*(rma + 4.0d0*pi*pre*rr**3)/ &
      &    (rr**2 - 2.0d0*rma*rr)
    f(3) =  4.0d0*pi*rr**2*rho0/(1.0d0-2.0d0*rma/rr)**0.5d0
    f(4) =  4.0d0*pi*rr**2*rho/(1.0d0-2.0d0*rma/rr)**0.5d0
    f(5) =  0.5d0*psi/rr*(1.0d0-1.0d0/(1.0d0-2.0d0*rma/rr)**0.5d0)
    f(6) =  4.0d0*pi*rr**2*rho/(1.0d0-2.0d0*rma/rr)**0.5d0/psi
  end if
!
  return
end subroutine test
!
! ---------------------------------------------------------
!
subroutine rk(neq, func, x0, xe, n, y0, yn, work)
!*********************************************************************
!    subroutine rk numerically integrates a system of neq
!    first order ordinary differential equations of the form
!            dy(i)/dx = f(x, y(1), ..., y(neq)),
!    by the classical runge-kutta formula.
!
!    parameters
!  === input ===
!    (1) neq: number of equations to be integrated
!    (2) func: subroutine func(x, y, f) to evaluate derivatives
!               f(i)=dy(i)/dx
!    (3) x0: initial value of independent variable
!    (4) xe: output point at which the solution is desired
!    (5) n: number of divisions
!       the interval (x0, xe) is divided into n subintervals
!       with the length (xe-x0)/n and in each subinterval
!       the classical runge-kutta formula is used.
!    (6) y0(i) (i=1, .., neq): initial value at x0
!  === output ===
!    (7) yn(i) (i=1, .., neq): approximate solution at xe
!  === other ===
!    (8) work(): two-dimentional array (size=(neq, 2)) to be
!                used inside rk
!    copyright: m. sugihara, november 15, 1989, v. 1
!*********************************************************************
  implicit none
  external func
  integer, intent(in)    :: neq, n
  real(8), intent(inout) :: x0, y0(neq), work(neq, 2)
  real(8), intent(in)    :: xe
  real(8), intent(out)   :: yn(neq)
  real(8) :: h
  integer :: i
  h = (xe - x0) / n
  do i = 1, n
    call rkstep(neq, func, x0, h, y0, yn, work(1, 1), work(1, 2))
    x0 = x0 + h
    y0(1:neq) = yn(1:neq)
  end do
  x0 = xe
  return
end subroutine rk
!
subroutine rkstep(neq, func, x, h, y0, yn, ak, w)
  implicit none
  real(8), parameter   :: a2 = 0.5d0, a3 = a2, &
                       &  b2 = 0.5d0, b3 = b2, &
                       &  c1 = 1.0d0 / 6, c2 = 1.0d0 / 3, c3 = c2, c4 = c1
  integer, intent(in)    :: neq
  real(8), intent(inout) :: x, h, y0(neq), ak(neq)
  real(8), intent(out)   :: yn(neq), w(neq)
  integer :: i
!
  call func(x, y0, ak)
  yn(1:neq) = y0(1:neq) + h * c1 * ak(1:neq)
  w(1:neq) = y0(1:neq) + h * b2 * ak(1:neq)
  call func(x + a2 * h, w, ak)
  do i = 1, neq
    yn(i) = yn(i) + h * c2 * ak(i)
  end do
  w(1:neq) = y0(1:neq) + h * b3 * ak(1:neq)
  call func(x + a3 * h, w, ak)
  do i = 1, neq
    yn(i) = yn(i) + h * c3 * ak(i)
  end do
  w(1:neq) = y0(1:neq) + h * ak(1:neq)
  call func(x + h, w, ak)
  do i = 1, neq
    yn(i) = yn(i) + h * c4 * ak(i)
  end do
  return
end subroutine rkstep
