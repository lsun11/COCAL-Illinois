subroutine fluidco(iter,iter_max,convf,cfvep,fmax0,emxemd,char)
!
  use phys_constant, only : pi
  use CB_fR_mesh_fluid, only : dr	!, drf, drfinv
  use def_metric_1D, only : alphf
  use CB_fR_param_physp, only : surr
  use def_matter_1D, only : emd, nnrg, emdc, ber, pinx, radi
  use grid_parameter_1D, only : nrf
  implicit none
!
  real(8) :: backf(0:nnrg)
  real(8) :: hut(0:nnrg), aloh(0:nnrg), hutp6(0:nnrg),psif4(0:nnrg)
  real(8) :: alpfc, emdfc, hhfc
  real(8), intent(inout) :: convf, cfvep, fmax0, emxemd
  integer, intent(inout) :: iter, iter_max
  integer :: iperr, ir, irerr,iterr
!
  character*3 char
!
! ----------------------------------------------------------------------
!
  dr = 1.0d0/dble(nrf)
!  drf = dr
!  drinv=1.0d0/dr
!  drfinv=drinv
!
! ----------------------------------------------------------------------
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!     Equations start. 
!
! ======================================================================
!
!     First integral of Euler equation.
!
! ======================================================================
!
! --  First integral of Euler equation.
!
  do ir = 0, nrf
    backf(ir) = emd(ir)
    emd(ir) = 1.0d0/(pinx+1.0d0) * (ber/alphf(ir) - 1.0d0)
  end do
!
  emd(0) = emdc
!
!  emdmx = 0.0d0
!  do ir = 0, nrf
!    emdmx = dmax1(emdmx,emd(ir))
!    if (emdmx == emd(ir)) nrfdmx = ir
!  end do
!
!  emd(0) = emdc
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
!     ### Iteration. ###
!     Set improved values for emden function and surface.
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!
! --- For parameters. 
!
  emd(nrf) = surr
  do ir = 0, nrf-1
    if (emd(ir) <= 0.0d0) emd(ir) = surr
  end do
!
  call flimpro(backf,emd,convf,emxemd,irerr,iterr,iperr,1)
  write(6,"(a22,es12.4,',  value =',es12.4,2x,3i4)") '  == emd   ==, error =', emxemd,&
     &              emd(irerr), irerr
!
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! --  Improving parameters and alpha and psi.  
!
  call paimproco_frg(ber,radi,convf,iter,fmax0,1)
!
! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!mark
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
end subroutine fluidco
