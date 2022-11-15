subroutine sourceterm_volume_int_bbh_2pot_test(sou)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi, alph
!  use def_matter, only : emdg
!  use def_matter_parameter, only : ber, radi, pinx
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint
  use make_array_3d
  implicit none
  real(long), pointer :: sou(:,:,:)
  real(8), pointer :: psi_grad2x(:,:,:), psi_grad2y(:,:,:), psi_grad2z(:,:,:)
  real(8), pointer :: alph_grad2x(:,:,:),alph_grad2y(:,:,:),alph_grad2z(:,:,:)
!  real(long) :: emdgc, rhogc, pregc, hhgc, utgc, rhoHc, zfac
  real(long) :: psigc, alpgc, alpsigc, psiinv
  real(long) :: dxpsi,dypsi,dzpsi,dxalph,dyalph,dzalph
  integer    :: irg, itg, ipg
!
  call alloc_array3d(psi_grad2x,1,nrg,1,ntg,1,npg)
  call alloc_array3d(psi_grad2y,1,nrg,1,ntg,1,npg)
  call alloc_array3d(psi_grad2z,1,nrg,1,ntg,1,npg)
  call alloc_array3d(alph_grad2x,1,nrg,1,ntg,1,npg)
  call alloc_array3d(alph_grad2y,1,nrg,1,ntg,1,npg)
  call alloc_array3d(alph_grad2z,1,nrg,1,ntg,1,npg)
!
  call grgrad_midpoint(psi,psi_grad2x,psi_grad2y,psi_grad2z)
  call grgrad_midpoint(alph,alph_grad2x,alph_grad2y,alph_grad2z)
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        psiinv = 1.0d0/psigc
        dxpsi  = psi_grad2x(irg,itg,ipg)
        dypsi  = psi_grad2y(irg,itg,ipg)
        dzpsi  = psi_grad2z(irg,itg,ipg)
        dxalph = alph_grad2x(irg,itg,ipg)
        dyalph = alph_grad2y(irg,itg,ipg)
        dzalph = alph_grad2z(irg,itg,ipg)
!
        sou(irg,itg,ipg) = -2.0d0*psiinv & 
        &                *(dxpsi*dxalph + dypsi*dyalph + dzpsi*dzalph)
      end do
    end do
  end do
!
  deallocate(psi_grad2x)
  deallocate(psi_grad2y)
  deallocate(psi_grad2z)
  deallocate(alph_grad2x)
  deallocate(alph_grad2y)
  deallocate(alph_grad2z)
end subroutine sourceterm_volume_int_bbh_2pot_test

