subroutine peos_sourcetermg(soupsg,souapg,soufpg)
!
  use phys_constant, only : nnrg, pi
  use def_metric_1D, only : psi, alph	!, alps
  use def_matter_1D, only : pinx, ber, radi, emdg
  use grid_parameter_1D, only : nrg, nrf
  implicit none
!
  real(8), intent(out) :: soupsg(0:nnrg), souapg(0:nnrg), soufpg(0:nnrg)
  real(8) :: x5(5), f5(5)
  real(8) :: alpgc, emdgc, hhgc, pregc, psigc, rhogc, rhoHc, rp2s, &
  &          utgc, zfac, ene
  integer :: irg
!
! --- Compute source terms for volume integrals.
!
!  intnum = 5
!  intp = 2
!  fac13 = 1.0d0/3.0d0
!  fac23 = 2.0d0/3.0d0
!  fac43 = 4.0d0/3.0d0
!  fac512 = 5.0d0/12.0d0
!
! --  GR coordinate
!
  soupsg(0:nrg)=0.0d0
  souapg(0:nrg)=0.0d0
  soufpg(0:nrg)=0.0d0
!
  do irg = 0, nrg
!
!    psiw = psi(irg)
!    alphw = alph(irg)
!    alpsw = alps(irg)
!    fphiw = fphi(irg)
!
! --  For fluid terms
!
    emdgc = emdg(irg)
    psigc =  psi(irg)
    alpgc = alph(irg)
!    rhogc = emdgc**pinx
!    pregc = rhogc*emdgc
!    hhgc  = 1.0d0 + (pinx+1.0d0)*emdgc
    call peos_q2hprho(emdgc, hhgc, pregc, rhogc, ene)
    utgc  = hhgc/ber
!
    zfac = 1.0d0
    if (emdgc <= 1.0d-15) zfac = 0.0d0
    rhoHc = hhgc*rhogc*(alpgc*utgc)**2 - pregc
    rp2s = 3.0d0*hhgc*rhogc*(alpgc*utgc)**2 - 2.0d0*hhgc*rhogc &
     &   + 5.0d0*pregc
!
! --  For psi.
    soupsg(irg) = - radi**2*2.0d0*pi*psigc**5*rhoHc*zfac
!
! --  For alpha*psi.
    souapg(irg) = + radi**2*2.0d0*pi*alpgc*psigc**5*rp2s*zfac
!
! --  For fphi.
    soufpg(irg) = + radi**2*2.0d0*pi*alpgc*psigc**5*rp2s*zfac
!
  end do
!
end subroutine peos_sourcetermg
