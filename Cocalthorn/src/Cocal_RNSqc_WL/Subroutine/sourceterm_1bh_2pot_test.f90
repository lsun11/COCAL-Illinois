subroutine sourceterm_1bh_2pot_test(sou)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use def_metric, only : psi, alph
!  use def_matter, only : emdg
!  use def_matter_parameter, only : ber, radi, pinx
  use interface_interpo_linear_type0
  use interface_grgrad_midpoint_r3rd_type0
  use interface_grgrad_midpoint_r4th_type0
  use make_array_3d
  implicit none
  real(long), pointer :: sou(:,:,:)
!  real(long) :: emdgc, rhogc, pregc, hhgc, utgc, rhoHc, zfac
  real(long) :: psigc, alpgc, alpsigc, psiinv
  real(long) :: dxpsi,dypsi,dzpsi,dxalph,dyalph,dzalph
  integer    :: irg, itg, ipg
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        call interpo_linear_type0(psigc,psi,irg,itg,ipg)
        call grgrad_midpoint_r3rd_type0(psi,dxpsi,dypsi,dzpsi,irg,itg,ipg,'bh')
        call grgrad_midpoint_r3rd_type0(alph,dxalph,dyalph,dzalph,&
                                                          &   irg,itg,ipg,'bh')
!       call grgrad_midpoint_r4th_type0(psi,dxpsi,dypsi,dzpsi,irg,itg,ipg,'bh')
!       call grgrad_midpoint_r4th_type0(alph,dxalph,dyalph,dzalph,&
!                                                         &   irg,itg,ipg,'bh')
        psiinv = 1.0d0/psigc
        sou(irg,itg,ipg) = -2.0d0*psiinv & 
        &                *(dxpsi*dxalph + dypsi*dyalph + dzpsi*dzalph)
      end do
    end do
  end do
!
end subroutine sourceterm_1bh_2pot_test
