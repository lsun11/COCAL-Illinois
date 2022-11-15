include '../../code/Module/phys_constant.f90'
include '../../code/Module/def_matter.f90'
include '../../code/Module/def_matter_parameter.f90'
include '../../code/Module/def_quantities.f90'
include '../../code/Module/grid_parameter.f90'
include '../../code/EOS/Module/def_teos_parameter.f90'
include '../../code/EOS/Module/def_peos_parameter.f90'
include '../../code/EOS/Subroutine/teos_initialize.f90'
include '../../code/EOS/Subroutine/teos_lookup.f90'
include '../../code/EOS/Subroutine/teos_h2qprho.f90'
include '../../code/EOS/Subroutine/teos_q2hprho.f90'
include '../../code/EOS/Subroutine/peos_lookup.f90'
include '../../code/EOS/Subroutine/peos_h2qprho.f90'
include '../../code/EOS/Subroutine/peos_q2hprho.f90'
include '../../code/EOS/Subroutine/peos2teos.f90'
include '../../code/EOS/Subroutine/stergioulas_rns2cocal.f90'
include './Module_TOV/def_TOV_quantities.f90'
include './Subroutine_TOV/rk.f90'
include './Subroutine_TOV/rkstep.f90'
include './Subroutine_TOV/equation_teos_isotr.f90'
include './Subroutine_TOV/printout_TOV_profile_teos_isotr.f90'
include './Subroutine_TOV/printout_TOV_quantities_teos_isotr.f90'
!
PROGRAM TOV
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
! --- Desctiption for parameters.
!
!    qini  : initial for emden function. 
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
!    yn_iso(1) : effective gravitational mass
!    yn_iso(2) : emden function
!    yn_iso(3) : proper rest mass 
!    yn_iso(4) : proper mass energy
!    yn_iso(5) : conformal factor (proper boundary condition not applied)
!    yn_iso(6) : ADM mass energy  (proper boundary condition not applied)
!    yn_iso(7) : lapse function
!    yn_iso(8) : Schwarzschild radial coordinate
!    adm  : ADM mass energy
!    compa : Compaction = YN(1)/XE = 
!              = Effective gravitational mass/
!                Radius of star in Schwartzscild radial coordinate 
!
! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
!
  use grid_parameter
  use def_matter_parameter, only : emdc, pinx
  use def_teos_parameter, only : emdini_gcm1
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
  real(8)  :: ar,xi,  q_xe, pre_xe, rho0_xe, ene_xe,  q1, h1, pre1, rho1, ene1
  real(8)  :: q1cgs, h1cgs, pre1cgs, rho1cgs, ene1cgs
!
  ermx      = 1.0d-10
  ermx_MorC = 1.0d-9
  itmx = 1000
!
  open (2, file='ovpar_teos.dat', status='old')
  open (9, file='ovaux.dat', status='unknown')
  read(2, '(2i5, 1es15.7)') nstep, ndiv, radiini
  read(2, '(2i5)') itype, iter_mode
  close(2)
!
! --  read_parameter in grid_parameter.f90
! --  set compactness and initial
!
  call read_parameter
!  call peos2teos(100)
  call stergioulas_rns2cocal
  call teos_initialize
  qini = emdini_gcm1
  call teos_q2hprho(qini, hini, pre, rho0, ene)
  write(6,'(a22,5e23.15)') "(qini,hini,P,rho,ene)=",qini, hini, pre, rho0, ene

!  call teos_h2qprho(hini, qini, pre, rho0, ene)
!  write(6,'(a22,5e23.15)') "(hini,qini,P,rho,ene)=",hini, qini, pre, rho0, ene
!stop

  dqdq = hini*0.1d0
!
! --- Iteration for compactness.
!
  erer = 1.0d0
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
      y0_iso(1) = 0.0d+0 ; y0_iso(2) = hini   ; y0_iso(3) = 0.0d+0
      y0_iso(4) = 0.0d+0 ; y0_iso(5) = 2.0d+0 ; y0_iso(6) = 0.0d+0
      y0_iso(7) = 0.0d+0 ; y0_iso(8) = 1.0d+0 ; y0_iso(9) = 0.0d+0
      y0_iso(10)= 0.0d+0 

      if (iter_radi > 1) then
        yr = yn_iso(9)/yn_iso(1)*(1.0d0-(1.0d0-2.0d0*yn_iso(1)/yn_iso(9))**0.5)
        ar = (1.0d0 - 2.0d0*yn_iso(1)/yn_iso(9))**0.5
        y0_iso(5) = yr*dexp(-yn_iso(7))
        y0_iso(8) = ar*dexp(-yn_iso(10))
        if ((radi_error <= ermx .and. erer < ermx_MorC) .or. (radi_error <= ermx .and. itype == 1)) then
          open (12, file='x0.las', status='unknown')
          call teos_h2qprho(hini, q1, pre1, rho1, ene1)
          rho1cgs = rho1*(c**6/((g**3)*solmas**2))
          pre1cgs = pre1*c**8/(g**3*solmas**2)
          ene1cgs = ene1*c**6/(g**3*solmas**2)
          write(12, '(17(es15.7))') 0.0d0, 0.0d0, hini, 0.0d0, 0.0d0, &
          &     y0_iso(5), 0.0d0, y0_iso(8), 0.0d0, rho1, rho1cgs, pre1, &
          &                        pre1cgs, ene1, ene1cgs, 0.0d0, 0.0d0

!          write(12,'(5es15.7)') x0, y0_iso(5), y0_iso(8), q1, hini
          close(12) 
        end if
      endif

      yn_iso(1) = 0.0d+0 ; yn_iso(2) = hini   ; yn_iso(3) = 0.0d+0 
      yn_iso(4) = 0.0d+0 ; yn_iso(5) = 1.0d+0 ; yn_iso(6) = 0.0d+0
      yn_iso(7) = 0.0d+0 ; yn_iso(8) = 1.0d+0 ; yn_iso(9) = 0.0d+0
      yn_iso(10)= 0.0d+0
      h = (xn - x0) / dble(nstep)
!
      do i = 1, nstep
        ii = i
        xe = x0 + h
        call rk(neq_iso, equation, x0, xe, ndiv, y0_iso, yn_iso, work_iso)
        if (yn_iso(2) <= 1.0d0) exit

        if ((radi_error <= ermx .and. itype == 1) .or. (radi_error <= ermx .and. erer < ermx_MorC)) then
          file_ctl = 'nm'
          if (i.eq.1)     file_ctl = 'op'
          if (i.eq.nstep) file_ctl = 'cl'
          call printout_TOV_profile_teos_isotr(file_ctl)
        end if
!        if (erer < ermx_MorC) then
!          if(i.eq.1)    open (11, file='rnsgra_1D_isotr.las', status='unknown')
!          call peos_h2qprho(yn_iso(2), q_xe, pre_xe, rho0_xe, ene_xe)
!          write(11,'(5es15.7)') xe/radi, yn_iso(5), yn_iso(8), q_xe, yn_iso(2)
!          if(i.eq.nstep)  close(11)
!        end if
      end do
!
      if (radi_error <= ermx .and. yn_iso(2)>=1.0d0) then
        if (itype == 0) exit
        if (itype == 1) then
          call printout_TOV_quantities_teos_isotr
          open(14, file='rnspar_add.dat', status='unknown')
          write(14, '(2es14.6, a15)') emdc, pinx, &
          & '  : emdc, pinx'
          write(14, '(2es14.6, a31)') restmass_sph, gravmass_sph, &
          & '  : restmass_sph, gravmass_sph'
          write(14, '(2es14.6, a50)') MoverR_sph, schwarz_radi_sph, &
          & '  : MoverR_sph,  schwarz_radi_sph  (G=c=M=1 unit)'
          write(6, '(/, /, 2es14.6, a15)') emdc, pinx, &
          & '  : emdc, pinx'
          write(6, '(2es14.6, a31)') restmass_sph, gravmass_sph, &
          & '  : restmass_sph, gravmass_sph'
          write(6, '(2es14.6, a50)') MoverR_sph, schwarz_radi_sph, &
          & '  : MoverR_sph,  schwarz_radi_sph  (G=c=M=1 unit)'
          close(14)
          call system("mv ovlas.dat ovlas.dat.tmp")
          call system("cat x0.las > ovlas.dat")
          call system("cat ovlas.dat.tmp >> ovlas.dat")
          call system("rm -f ovlas.dat.tmp")
          stop
        end if
      end if
!
! ---- hunting for precise radius.
!
      if (yn_iso(2) <= 1.0d0) then
        if (iter_radi == 1) write(6, *) ' bad initial radius '
        radi = radi - dr
        dr_back = dr
        write(9, *)' back ', iter_radi, ii, radi, radi_error, yn_iso(2)
        call system("rm -f ./ovlas.dat")
      else
        radi = radi + 0.2*dr_back
        radi_error = dr_back/radi
        dr = 0.2d0*dr_back
        write(9, *)' hunt ', iter_radi, ii, radi, radi_error, yn_iso(2)
      end if
      if (iter_radi == itmx)     stop ' iter_radi '
!
! ---- hunting end.
!
    end do
!
! ---- Searching for a certain value of compactness.
!
    yr = yn_iso(9)/yn_iso(1)*(1.0d0-(1.0d0-2.0d0*yn_iso(1)/yn_iso(9))**0.5)
!    write(6,'(a59,es15.7)')  "psi(R) from exact formula in Schwarzschild coordinates yr: ", yr
!    write(6,'(a59,es15.7)')  "psi(R) calculated, yn_iso(5)                             : ", yn_iso(5)
 
    adm = yn_iso(5)*yn_iso(6)/yr
    compa = yn_iso(1)/yn_iso(9)
!
    if (iter_mode.eq.0) then
      Q_current = compa
      Q_target  = MoverR_sph
    else if (iter_mode.eq.1) then
      Q_current = yn_iso(3)
      Q_target  = restmass_sph
    else if (iter_mode.eq.2) then
      Q_current = yn_iso(1)
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
    if (erer < ermx_MorC) then
      call printout_TOV_quantities_teos_isotr
      open(14, file='rnspar_add.dat', status='unknown') 
      write(14, '(2es14.6, a15)') emdc, pinx, &
        & '  : emdc, pinx'
      write(14, '(2es14.6, a31)') restmass_sph, gravmass_sph, &
        & '  : restmass_sph, gravmass_sph'
      write(14, '(2es14.6, a50)') MoverR_sph, schwarz_radi_sph, &
        & '  : MoverR_sph,  schwarz_radi_sph  (G=c=M=1 unit)'
      write(6, '(/, /, 2es14.6, a15)') emdc, pinx, &
        & '  : emdc, pinx'
      write(6, '(2es14.6, a31)') restmass_sph, gravmass_sph, &
        & '  : restmass_sph, gravmass_sph'
      write(6, '(2es14.6, a50)') MoverR_sph, schwarz_radi_sph, &
        & '  : MoverR_sph,  schwarz_radi_sph  (G=c=M=1 unit)'
      close(14)
      call system("mv ovlas.dat ovlas.dat.tmp")
      call system("cat x0.las > ovlas.dat")
      call system("cat ovlas.dat.tmp >> ovlas.dat")
      call system("rm -f ovlas.dat.tmp")
      stop
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
  end do
!
  stop
end PROGRAM TOV
!
function fn_lagint(x,y,v)
  implicit none
  real(8), intent(in) :: x(4),y(4), v
  real(8) :: dx12, dx13, dx14, dx21, dx23, dx24, dx31, dx32, dx34, &
  &          dx41, dx42, dx43, wex1, wex2, wex3, wex4, &
  &          xv1, xv2, xv3, xv4, fn_lagint
!
  dx12 = x(1) - x(2)
  dx13 = x(1) - x(3)
  dx14 = x(1) - x(4)
  dx23 = x(2) - x(3)
  dx24 = x(2) - x(4)
  dx34 = x(3) - x(4)
  dx21 = - dx12
  dx31 = - dx13
  dx32 = - dx23
  dx41 = - dx14
  dx42 = - dx24
  dx43 = - dx34
  xv1 = v - x(1)
  xv2 = v - x(2)
  xv3 = v - x(3)
  xv4 = v - x(4)
  wex1 = xv2*xv3*xv4/(dx12*dx13*dx14)
  wex2 = xv1*xv3*xv4/(dx21*dx23*dx24)
  wex3 = xv1*xv2*xv4/(dx31*dx32*dx34)
  wex4 = xv1*xv2*xv3/(dx41*dx42*dx43)
!
  fn_lagint = wex1*y(1) + wex2*y(2) + wex3*y(3) + wex4*y(4)
!
end function fn_lagint

!
function fn_lagint_2nd(x,y,v)
  implicit none
  real(8), intent(in) :: x(2),y(2), v
  real(8) :: dx12, dx21, wex1, wex2, xv1, xv2, fn_lagint_2nd
!
  dx12 = x(1) - x(2)
  dx21 = - dx12
  xv1 = v - x(1)
  xv2 = v - x(2)
  wex1 = xv2/dx12
  wex2 = xv1/dx21
!
  fn_lagint_2nd = wex1*y(1) + wex2*y(2)
!
end function fn_lagint_2nd

