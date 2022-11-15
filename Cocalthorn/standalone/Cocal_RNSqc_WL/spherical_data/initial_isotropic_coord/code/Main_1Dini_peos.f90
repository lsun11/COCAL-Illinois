!______________________________________________
include 'include_modulefiles_1D_peos.f90'
include 'include_subroutines_1D_peos.f90'
!______________________________________________
!
!              Main Program
!______________________________________________
PROGRAM GR_ISOTROPIC_NS_1D
!______________________________________________
! Spherical neutron star in GR 
! calculated in isotropic coordinate using KEH iteration.
!______________________________________________
!
  use phys_constant, only : nnrg
  use def_metric_1D, only : alphf, psif, psi, alps, alph
  use CB_fR_param_physp, only : surr
  use CB_fR_param_iomis, only : itype, numseq
  use grid_parameter_1D, only : iter_max, rgmid, nrgin, eps, emdc, nrg, &
  &                          conv_ini, conv0_gra, conv0_den, &
  &                          chgra, chope, chrot
  use def_matter_1D, only : ber, emd, emdg, radi
  use coordinate_grav_r_1D, only : rg, hrg
!
  implicit none
!
!ref  common / physp / surr
!ref  common / flmas / ahores
!ref  common / iomis / numseq, itype
!
  real(8) :: soupsg(0:nnrg), souapg(0:nnrg), soufpg(0:nnrg)
  real(8) :: potg(0:nnrg), backg(0:nnrg)
  real(8) :: cfvep, convf, dsouapou, dsoufpou, dsoupsou, emxemd, &
  &          epsmax, fcon, fffac, fmax0, &
  &          epsmaxg, souapou, soufpou, soupsou
  real(8) :: emdgc, hhgc, pregc, rhogc, ene
  integer :: ii, irger, iseq, istat, iter, ir
!
  character*4 char
  character*1 chbhs
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  open(26,file='frgaux2.dat',status='unknown')
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
  write(6,*)' coordinates etc. '
!
  call coordinate_patch_kit_grav_1D
!
  itype = 0
  chbhs = 'c'
  numseq = 1
!
! +++ Coordinates
!
!###  Fluid coordinate
!
  char = chrot//chbhs//chgra//chope
!
! ----------------------------------------------------------------------
! --- Initial.
! ----------------------------------------------------------------------
!
  istat = itype
  iseq  = 0
  call subio_1D(istat,iseq,char)
  call subio_fluid_1D(istat,iseq,char)
!
  call peos_initialize
! ----------------------------------------------------------------------
! --- Calcuration of a sequence.
! ----------------------------------------------------------------------
!
! --  Initial might be before nonv**.
!
  do iseq = 1, numseq
!
    surr = 10.0d0**(idnint(dlog10(1.0d-04*emdc)))
!!    surr = emdc**1.0d-20
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Iteration starts.  
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
    do iter = 1, iter_max
!
      epsmax = 0.0d0
      fmax0 = 0.0d0
      fcon = conv_ini*dble(iter)
      fffac = dmin1(conv0_gra, fcon)
      convf = dmin1(conv0_den, fcon)
!
! --- Start computation.  
!
      write(6,*) ' #### Iteration NO. = ', iter
!
! ==  Interpolations from fluid.
!
      call al2alps
      call fl2gr(emd,emdg)
!
! --- Compute source terms.
!
      call asymptopia_alps(soupsou,souapou,dsoupsou,dsouapou)
      call peos_sourcetermg(soupsg,souapg,soufpg)
!
! --- Call Poisson solver to compute Green's formula.
!
      do ii = 1, 2
!     do ii = 1, 3
!
        potg(0:nrg) = 0.0d0
        if (ii == 1) potg(0:nrg) = psi(0:nrg)
        if (ii == 2) potg(0:nrg) = alps(0:nrg)
!       if (ii == 3) potg(0:nrg) = fphi(0:nrg)
        backg(0:nrg) = potg(0:nrg)
!
!
        if (ii == 1) call poisol(soupsou,dsoupsou,soupsg,potg,0,0,0,ii)
        if (ii == 2) call poisol(souapou,dsouapou,souapg,potg,0,0,0,ii)
!       if (ii == 3) call poisol(soufpou,dsoufpou,soufpg,potg,0,0,0,ii)
!
        if (ii == 1) write(6,*)'--- max error in GR coordinate ---'
!
        if (ii == 1) write(6,*) '=== psi     === '
        if (ii == 2) write(6,*) '=== alpha   === '
!       if (ii == 3) write(6,*) '=== fphi    === '
        call frgerror(backg,potg,epsmaxg,irger)
        epsmax = dmax1(epsmaxg,epsmax)
        write(6,"(' GR error ',i5,3es12.4)") &
         & irger, backg(irger), potg(irger), epsmaxg
!
! --  Update variables.
!
        call potimpro(potg,backg,fffac)
!
        if (ii == 1) psi(0:nrg) = potg(0:nrg)
        if (ii == 2) alps(0:nrg) = potg(0:nrg)
!       if (ii == 3) fphi(0:nrg) = potg(0:nrg)
!
        if (ii == 2) call alps2al
        call peos_paimproco_frg(ber,radi,convf,iter,fmax0,0)
        epsmax = dmax1(fmax0,epsmax)
!
      end do
!
      call alps2al
      call gr2fl(alph,alphf)
      call gr2fl(psi,psif)
!
! --  Solve First integral of Euler eq.
!
      call peos_fluidco(iter,iter_max,convf,cfvep,fmax0,emxemd,char)
      epsmax = dmax1(fmax0,emxemd,epsmax)
!
      call peos_restmass
      call peos_admmass
!
      write(6,'(a22,2es12.4)')'  -- epsmax --, error =', epsmax
!
! --- Untill convergence is made, go to 2444.
      if (epsmax <= eps) exit
    end do
!
! --- Print out the intermediate state.
    if (epsmax > eps) then
      istat = 2
    else
      istat = 1
    end if
!
!    open(26,file='frgaux2.dat',status='unknown')
!    call peos_sourcetermg(soupsg,souapg,soufpg)
!    do ir = 0, nrg
!      write(26,*) rg(ir), soupsg(ir)
!!      emdgc = emdg(ir)
!!      call peos_q2hprho(emdgc, hhgc, pregc, rhogc, ene)
!!      write(26,*) rg(ir), emdgc, hhgc, pregc, rhogc, ene
!    end do
!      write(26,*) ' '
!!    call halfsou(soupsg)
!    call interpo_fl2grmidpoint_1D(soupsg)
!!    do ir = 0, nrg
!    do ir = 1, nrg
!      write(26,*) hrg(ir), soupsg(ir)
!!      emdgc = emdg(ir)
!!      call peos_q2hprho(emdgc, hhgc, pregc, rhogc, ene)
!!      write(26,*) rg(ir), emdgc, hhgc, pregc, rhogc, ene
!    end do
!    close(26)
!!
! --- Print out the converged state.
    call subio_1D(istat,iseq,char)
    call subio_fluid_1D(istat,iseq,char)
    call peos_physq_fR(istat,iseq)
    call printout_plot
!
    if (epsmax > eps) then
      write(6,*) ' **iter** '
      write(6,*) &
       & ' ## Iteration did not converge - repeat from .nxt data ## '
    else
      write(6,*) ' ## Converged ## '
    end if
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --- Calcuration of deformation sequence.
!
    write(6,*)' Sequence No.= ',iseq
!
  end do
!
! --- end of the program
!
END PROGRAM GR_ISOTROPIC_NS_1D
!______________________________________________
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
!
