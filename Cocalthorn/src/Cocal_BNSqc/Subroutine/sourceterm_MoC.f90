subroutine sourceterm_MoC(souvec,sou)
  use phys_constant, only : long, pi
  use grid_parameter, only : nrg, ntg, npg
  use coordinate_grav_r, only : hrg, rg
  use trigonometry_grav_theta, only : hsinthg,  hcosthg
  use trigonometry_grav_phi,   only : hsinphig, hcosphig
  use def_metric, only : psi, alph, tfkijkij, bvxd, bvyd, bvzd, tfkij
  use def_matter, only : emdg
  use def_matter_parameter, only : ome, ber, radi, pinx
  use def_vector_x,   only : hvec_xg
  use def_vector_phi, only : hvec_phig
  use make_array_3d
  use interface_grgrad_midpoint
  use interface_interpo_linear_type0
  implicit none
  real(long), pointer :: sou(:,:,:) 
  real(long), pointer :: souvec(:,:,:,:)
  real(long), pointer :: fnc2(:,:,:)
  real(long), pointer :: grad2x(:,:,:), grad2y(:,:,:), grad2z(:,:,:)
  real(long) :: vphig(3), xxx, yyy, zzz, san, san2
  real(long) :: emdgc, rhogc, pregc, hhgc, utgc, oterm, zfac, rjj
  real(long) :: psigc, alpgc, fnc2gc, afnc2inv, dxafn, dyafn, dzafn
  real(long) :: bvxdgc, bvydgc, bvzdgc, tfkax, tfkay, tfkaz
  integer :: ii, irg, itg, ipg
!
  call alloc_array3d(fnc2,0,nrg,0,ntg,0,npg)
  call alloc_array3d(grad2x,1,nrg,1,ntg,1,npg)
  call alloc_array3d(grad2y,1,nrg,1,ntg,1,npg)
  call alloc_array3d(grad2z,1,nrg,1,ntg,1,npg)
!
  san = 1.0d0/3.0d0
  san2= 2.0d0/3.0d0
!
! --- Source terms of Momentum constraint 
! --- for computing shift.  
!
  do ipg = 0, npg
    do itg = 0, ntg
      do irg = 0, nrg
        fnc2(irg,itg,ipg) = psi(irg,itg,ipg)**6/alph(irg,itg,ipg)
      end do
    end do
  end do
!
  call grgrad_midpoint(fnc2,grad2x,grad2y,grad2z)
!
  do ii = 1, 3
    do ipg = 1, npg
      do itg = 1, ntg
        do irg = 1, nrg
          call interpo_linear_type0(emdgc,emdg,irg,itg,ipg)
          call interpo_linear_type0(psigc,psi,irg,itg,ipg)
          call interpo_linear_type0(alpgc,alph,irg,itg,ipg)
          call interpo_linear_type0(fnc2gc,fnc2,irg,itg,ipg)
          call interpo_linear_type0(bvxdgc,bvxd,irg,itg,ipg)
          call interpo_linear_type0(bvydgc,bvyd,irg,itg,ipg)
          call interpo_linear_type0(bvzdgc,bvzd,irg,itg,ipg)
          afnc2inv = alpgc/fnc2gc
          dxafn = grad2x(irg,itg,ipg)
          dyafn = grad2y(irg,itg,ipg)
          dzafn = grad2z(irg,itg,ipg)
          tfkax  = tfkij(irg,itg,ipg,ii,1)
          tfkay  = tfkij(irg,itg,ipg,ii,2)
          tfkaz  = tfkij(irg,itg,ipg,ii,3)
!
          zfac = 1.0d0
          if (emdgc <= 1.0d-15) then
            emdgc = 1.0d-15
            zfac  = 0.0d0
          end if
          rhogc = emdgc**pinx
          pregc = rhogc*emdgc
          hhgc  = 1.0d0 + (pinx+1.0d0)*emdgc
          utgc  = hhgc/ber
          vphig(1) = hvec_phig(irg,itg,ipg,1)
          vphig(2) = hvec_phig(irg,itg,ipg,2)
          vphig(3) = hvec_phig(irg,itg,ipg,3)
          oterm = 0.0d0
          if (ii == 1) oterm = bvxdgc + ome*vphig(1)
          if (ii == 2) oterm = bvydgc + ome*vphig(2)
          if (ii == 3) oterm = bvzdgc + ome*vphig(3)
          rjj = hhgc*rhogc*alpgc*utgc**2*psigc**4*oterm
!
          souvec(irg,itg,ipg,ii) = - 2.0d0*afnc2inv &
        &   *(tfkax*dxafn + tfkay*dyafn + tfkaz*dzafn) &
        &   + radi**2*16.0d0*pi*alpgc*rjj*zfac
        end do
      end do
    end do
  end do     
!
! --- For a function B.
!
  do ipg = 1, npg
    do itg = 1, ntg
      do irg = 1, nrg
        xxx = hvec_xg(irg,itg,ipg,1)
        yyy = hvec_xg(irg,itg,ipg,2)
        zzz = hvec_xg(irg,itg,ipg,3)
        sou(irg,itg,ipg) = xxx*souvec(irg,itg,ipg,1)  &
      &                  + yyy*souvec(irg,itg,ipg,2)  &
      &                  + zzz*souvec(irg,itg,ipg,3)
      end do
    end do
  end do
!
  deallocate(fnc2)
  deallocate(grad2x)
  deallocate(grad2y)
  deallocate(grad2z)
end subroutine sourceterm_MoC
